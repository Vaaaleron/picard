/*
 * The MIT License
 *
 * Copyright (c) 2015 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package picard.analysis;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFileWalker;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.SequenceUtil;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.Semaphore;
import java.util.concurrent.atomic.AtomicBoolean;

/**
 * Super class that is designed to provide some consistent structure between
 * subclasses that simply iterate once over a coordinate sorted BAM and collect
 * information from the records as the go in order to produce some kind of
 * output.
 *
 * @author Tim Fennell
 */
public abstract class SinglePassSamProgram extends CommandLineProgram {
	@Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "Input SAM or BAM file.")
	public File INPUT;

	@Option(shortName = "O", doc = "File to write the output to.")
	public File OUTPUT;

	@Option(doc = "If true (default), then the sort order in the header file will be ignored.", shortName = StandardOptionDefinitions.ASSUME_SORTED_SHORT_NAME)
	public boolean ASSUME_SORTED = true;

	@Option(doc = "Stop after processing N reads, mainly for debugging.")
	public long STOP_AFTER = 0;

	private static final Log log = Log.getInstance(SinglePassSamProgram.class);

	/**
	 * Final implementation of doWork() that checks and loads the input and
	 * optionally reference sequence files and the runs the sublcass through the
	 * setup() acceptRead() and finish() steps.
	 */
	@Override
	protected final int doWork() {
		makeItSo(INPUT, REFERENCE_SEQUENCE, ASSUME_SORTED, STOP_AFTER, Arrays.asList(this));
		return 0;
	}

	public static void makeItSo(final File input, final File referenceSequence, final boolean assumeSorted,
			final long stopAfter, final Collection<SinglePassSamProgram> programs) {

		// Setup the standard inputs
		IOUtil.assertFileIsReadable(input);
		final SamReader in = SamReaderFactory.makeDefault().referenceSequence(referenceSequence).open(input);

		// Optionally load up the reference sequence and double check sequence
		// dictionaries
		final ReferenceSequenceFileWalker walker;
		if (referenceSequence == null) {
			walker = null;
		} else {
			IOUtil.assertFileIsReadable(referenceSequence);
			walker = new ReferenceSequenceFileWalker(referenceSequence);

			if (!in.getFileHeader().getSequenceDictionary().isEmpty()) {
				SequenceUtil.assertSequenceDictionariesEqual(in.getFileHeader().getSequenceDictionary(),
						walker.getSequenceDictionary());
			}
		}

		// Check on the sort order of the BAM file
		{
			final SortOrder sort = in.getFileHeader().getSortOrder();
			if (sort != SortOrder.coordinate) {
				if (assumeSorted) {
					log.warn("File reports sort order '" + sort + "', assuming it's coordinate sorted anyway.");
				} else {
					throw new PicardException("File " + input.getAbsolutePath() + " should be coordinate sorted but "
							+ "the header says the sort order is " + sort + ". If you believe the file "
							+ "to be coordinate sorted you may pass ASSUME_SORTED=true");
				}
			}
		}

		// Call the abstract setup method!
		boolean anyUseNoRefReads = false;
		for (final SinglePassSamProgram program : programs) {
			program.setup(in.getFileHeader(), input);
			anyUseNoRefReads = anyUseNoRefReads || program.usesNoRefReads();
		}

		final ProgressLogger progress = new ProgressLogger(log);

		ExecutorService service = Executors.newFixedThreadPool(2);

		Semaphore sem = new Semaphore(1);

		final int QUEUE_CAPACITY = 10;
		
		Object monitor = new Object();

		class Worker implements Runnable {

			final List<Object[]> poisonPill = Collections.emptyList();

			BlockingQueue<List<Object[]>> queue = new LinkedBlockingQueue<List<Object[]>>(QUEUE_CAPACITY);

			AtomicBoolean working = new AtomicBoolean(false);

			@Override
			public void run() {
				while (true) {
					try {
						final List<Object[]> tmpPairs = queue.take();
						
						if (tmpPairs.isEmpty()) {
							return;
						}
						
						sem.acquire();
						service.submit(new Runnable() {

							@Override
							public void run() {
								synchronized (monitor) {
									working.set(true);
									monitor.notify();
								}
								for (Object[] objects : tmpPairs) {
									final SAMRecord rec = (SAMRecord) objects[0];
									final ReferenceSequence ref = (ReferenceSequence) objects[1];

									for (final SinglePassSamProgram program : programs) {
										program.acceptRead(rec, ref);
									}

								}

								synchronized (monitor) {
									working.set(false);
									monitor.notify();
								}
								sem.release();
							}
						});

					} catch (InterruptedException e) {
						e.printStackTrace();
					}
				}
			}

			public void submitData(List<Object[]> data) {
				try {
					queue.put(data);
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}

			public void stop() {
				try {
					queue.put(poisonPill);
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
		}

		// Starting ExecutorService
		Worker worker = new Worker();
		service.execute(worker);

		Iterator<SAMRecord> it = in.iterator();

		// finding optimal MAX_PAIRS value
		// first test
		SAMRecord rec = it.next();

		ReferenceSequence ref;
		if (walker == null || rec.getReferenceIndex() == SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX) {
			ref = null;
		} else {
			ref = walker.get(rec.getReferenceIndex());
		}

		int MAX_PAIRS = 1;
		List<Object[]> pairs = new ArrayList<>(MAX_PAIRS);

		pairs.add(new Object[] { rec, ref });

		long startSub = System.nanoTime();

		worker.submitData(pairs);

		synchronized (monitor) {
			while (!worker.working.get()) {
				try {
					monitor.wait();
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
		}
		
		long startWrk = System.nanoTime();

		synchronized (monitor) {
			while (worker.working.get()) {
				try {
					monitor.wait();
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
		}
		long stopWrk = System.nanoTime();

		progress.record(rec);

		long submitting1 = startWrk - startSub;
		long working1 = stopWrk - startWrk;

		// finding optimal MAX_PAIRS value
		// second test
		if (((stopAfter > 2) || (stopAfter == 0)) && (submitting1 > working1)) {
			MAX_PAIRS = 2;
			pairs = new ArrayList<>(MAX_PAIRS);

			while (it.hasNext() && pairs.size() < MAX_PAIRS) {
				rec = it.next();

				if (walker == null || rec.getReferenceIndex() == SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX) {
					ref = null;
				} else {
					ref = walker.get(rec.getReferenceIndex());
				}

				progress.record(rec);
				pairs.add(new Object[] { rec, ref });

				if (!anyUseNoRefReads && rec.getReferenceIndex() == SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX) {
					break;
				}
			}

			startSub = System.nanoTime();

			worker.submitData(pairs);

			synchronized (monitor) {
				while (!worker.working.get()) {
					try {
						monitor.wait();
					} catch (InterruptedException e) {
						e.printStackTrace();
					}
				}
			}
			startWrk = System.nanoTime();

			synchronized (monitor) {
				while (worker.working.get()) {
					try {
						monitor.wait();
					} catch (InterruptedException e) {
						e.printStackTrace();
					}
				}
			}
			stopWrk = System.nanoTime();


			long submitting2 = startWrk - startSub;
			long working2 = stopWrk - startWrk;

			long wrkOneElement = working2 - working1;
			long subOneElement = submitting2 - submitting1;
			long subQueue = submitting1 - subOneElement;

			if (wrkOneElement > subOneElement) {
				MAX_PAIRS = (int) (subQueue / (wrkOneElement - subOneElement));
			} else {
				MAX_PAIRS = 10000;
			}

		}

		pairs = new ArrayList<>(MAX_PAIRS);

		while (it.hasNext()) {
			// See if we need to terminate early?
			if (stopAfter > 0 && progress.getCount() >= stopAfter) {
				break;
			}

			// And see if we're into the unmapped reads at the end
			if (!anyUseNoRefReads && rec.getReferenceIndex() == SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX) {
				break;
			}

			rec = it.next();

			if (walker == null || rec.getReferenceIndex() == SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX) {
				ref = null;
			} else {
				ref = walker.get(rec.getReferenceIndex());
			}

			progress.record(rec);
			pairs.add(new Object[] { rec, ref });
			if (pairs.size() < MAX_PAIRS) {
				continue;
			}

			worker.submitData(pairs);

			pairs = new ArrayList<>(MAX_PAIRS);
		}
		
		if (pairs.size() > 0) {
			worker.submitData(pairs);
		}

		service.shutdown();
		worker.stop();

		CloserUtil.close(in);

		for (final SinglePassSamProgram program : programs) {
			program.finish();
		}
	}

	/**
	 * Can be overriden and set to false if the section of unmapped reads at the
	 * end of the file isn't needed.
	 */
	protected boolean usesNoRefReads() {
		return true;
	}

	/**
	 * Should be implemented by subclasses to do one-time initialization work.
	 */
	protected abstract void setup(final SAMFileHeader header, final File samFile);

	/**
	 * Should be implemented by subclasses to accept SAMRecords one at a time.
	 * If the read has a reference sequence and a reference sequence file was
	 * supplied to the program it will be passed as 'ref'. Otherwise 'ref' may
	 * be null.
	 */
	protected abstract void acceptRead(final SAMRecord rec, final ReferenceSequence ref);

	/** Should be implemented by subclasses to do one-time finalization work. */
	protected abstract void finish();

}

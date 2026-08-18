/*
 * 
 * “Commons Clause” License Condition v1.0
 * 
 * The Software is provided to you by the Licensor under the License, 
 * as defined below, subject to the following condition.
 * 
 * Without limiting other conditions in the License, the grant of rights under the License 
 * will not include, and the License does not grant to you, the right to Sell the Software.
 * 
 * For purposes of the foregoing, “Sell” means practicing any or all of the rights granted 
 * to you under the License to provide to third parties, for a fee or other consideration 
 * (including without limitation fees for hosting or consulting/ support services related to 
 * the Software), a product or service whose value derives, entirely or substantially, from the 
 * functionality of the Software. Any license notice or attribution required by the License 
 * must also include this Commons Clause License Condition notice.
 * 
 * Software: genestrip
 * 
 * License: Apache 2.0
 * 
 * Licensor: Daniel Pfeifer (daniel.pfeifer@progotec.de)
 * 
 */
package org.metagene.genestrip.goals.refseq;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.atomic.AtomicInteger;

import org.junit.Test;
import org.metagene.genestrip.DefaultExecutionContext;
import org.metagene.genestrip.ExecutionContext;
import org.metagene.genestrip.GSCommon;
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.refseq.AbstractRefSeqFastaReader;
import org.metagene.genestrip.tax.TaxTree;

/**
 * Tests that a reader goal can read a second time in the same process.
 * <p>
 * It could not, and the way it failed was quiet. A consumer thread loops on the queue until the flag
 * {@code dump} tells it to stop, and a pass used to end by setting that flag: every subclass did so
 * itself, from the {@code finally} that released its own state. Nothing ever put it back, so a goal
 * that had read once was finished for good -- asked again, it queued every file, its consumers saw
 * the flag still set and stopped before touching one, and the pass returned in seconds having read
 * nothing at all. No exception, no warning, just an empty result that looks like an empty input.
 * <p>
 * A pass now ends its own consumers with a sentinel per consumer and leaves the flag alone, so the
 * flag means only what its name says. This goal therefore does what the reworked subclasses do:
 * it reads and keeps its result, and ends nothing by hand.
 * <p>
 * That is what a linkage sweep ran into on its second linkage, twice, and it is why the invariant is
 * pinned here rather than left to the goals that happen to depend on it.
 */
public class FastaReaderGoalSecondPassTest {
    private static final TaxTree.TaxIdNode NODE = new TaxTree.TaxIdNode("9606");

    /** A reader that only counts the regions it is handed, so a pass is cheap. */
    private static final class CountingReader extends AbstractRefSeqFastaReader {
        private final AtomicInteger regions;

        CountingReader(Set<TaxTree.TaxIdNode> taxNodes, AtomicInteger regions,
                       StringLong2DigitTrie regionsPerTaxid) {
            super(4096, taxNodes, null, 31, Integer.MAX_VALUE, null, Long.MAX_VALUE, 1, false,
                    regionsPerTaxid);
            this.regions = regions;
        }

        @Override
        protected void infoLine() {
            super.infoLine();
            regions.incrementAndGet();
        }

        @Override
        protected void dataLine() {
        }
    }

    /** A goal that reads the FASTAs it is given and ends its pass the way the real ones do. */
    private static final class CountingGoal extends FastaReaderGoal<Integer, GSProject> {
        private final AtomicInteger regions = new AtomicInteger();

        @SuppressWarnings("unchecked")
        CountingGoal(GSProject project, ExecutionContext bundle,
                     ObjectGoal<Map<File, TaxTree.TaxIdNode>, GSProject> additionalGoal) {
            super(project, GSGoalKey.FILLSIZE, bundle, constant(project, GSGoalKey.CATEGORIES,
                            Collections.<org.metagene.genestrip.refseq.RefSeqCategory>emptySet()),
                    constant(project, GSGoalKey.TAXNODES, Collections.<TaxTree.TaxIdNode>emptySet()),
                    null, additionalGoal, false);
        }

        @Override
        protected AbstractRefSeqFastaReader createFastaReader(
                AbstractRefSeqFastaReader.StringLong2DigitTrie regionsPerTaxid) {
            return new CountingReader(Collections.singleton(NODE), regions, regionsPerTaxid);
        }

        @Override
        protected void doMakeThis() {
            try {
                readFastas();
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
            // Nothing else to do: the pass has ended its own consumers. Subclasses used to call
            // cleanUpThreads() here, which is what left the goal unusable afterwards.
            set(regions.get());
        }

        int readAgain() {
            regions.set(0);
            doMakeThis();
            return regions.get();
        }
    }

    /** An object goal that simply holds a value. */
    private static <T> ObjectGoal<T, GSProject> constant(GSProject project, GSGoalKey key, T value) {
        return new ObjectGoal<T, GSProject>(project, key) {
            @Override
            protected void doMakeThis() {
                set(value);
            }
        };
    }

    @Test
    public void aGoalThatHasReadOnceCanReadAgain() throws IOException {
        File dir = Files.createTempDirectory("fastareader").toFile();
        dir.deleteOnExit();
        Map<File, TaxTree.TaxIdNode> fastas = new LinkedHashMap<>();
        for (int i = 0; i < 4; i++) {
            File fasta = new File(dir, "g" + i + ".fna");
            Files.write(fasta.toPath(), (">g" + i + " test\nACGTACGTACGTACGTACGTACGTACGTACGTACGT\n")
                    .getBytes(StandardCharsets.UTF_8));
            // A node with the file keeps the reader off the accession map, which this test
            // has no reason to build: what is under test is the pass, not the mapping.
            fastas.put(fasta, NODE);
        }

        GSProject project = new GSProject(new GSCommon(dir), "test", null, null, null, null, null,
                null, null, null, null, false);
        project.initConfigParam(GSConfigKey.PROGRESS_BAR, false);
        // No main thread handed over: ending a pass interrupts it along with the consumers, and this
        // test drives two passes from the thread that would then be carrying the interrupt into the
        // second one's blockingQueue.put().
        ExecutionContext bundle = new DefaultExecutionContext(null, 2, 1000000);

        CountingGoal goal = new CountingGoal(project, bundle, constant(project, GSGoalKey.ADD_FASTAS, fastas));
        goal.make();
        assertEquals("the first pass must read every file", 4, goal.get().intValue());
        // The second pass is the one that used to read nothing while reporting success.
        assertEquals("a second pass must read every file again", 4, goal.readAgain());
    }
}

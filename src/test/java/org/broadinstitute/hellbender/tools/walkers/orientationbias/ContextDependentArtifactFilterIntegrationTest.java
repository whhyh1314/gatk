package org.broadinstitute.hellbender.tools.walkers.orientationbias;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.testng.annotations.Test;

import java.io.File;

import static org.testng.Assert.*;

/**
 * Created by tsato on 8/1/17.
 */
public class ContextDependentArtifactFilterIntegrationTest extends CommandLineProgramTest {

    @Test
    public void test() {
        final File output = createTempFile("contamination", ".table");
        final File mutectVcf = new File(DREAM_VCFS_DIR, "sample_1.vcf");
        final File dream1Bam = new File(DREAM_BAMS_DIR, "tumor_1.bam");
        final String[] args = {
                "-R", b37_reference_20_21,
                "-I", dream1Bam.getAbsolutePath(),
                "-O", output.getAbsolutePath()
        };


        runCommandLine(args);

    }

}
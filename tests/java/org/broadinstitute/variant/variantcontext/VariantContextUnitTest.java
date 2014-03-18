/*
* Copyright (c) 2012 The Broad Institute
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.variant.variantcontext;


// the imports for unit testing.

import org.broadinstitute.variant.VariantBaseTest;
import org.testng.annotations.BeforeSuite;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import org.testng.Assert;

import java.lang.reflect.Array;
import java.util.*;


public class VariantContextUnitTest extends VariantBaseTest {
    Allele A, Aref, C, T, Tref;
    Allele del, delRef, ATC, ATCref;

    // A [ref] / T at 10
    String snpLoc = "chr1";
    int snpLocStart = 10;
    int snpLocStop = 10;

    // - / ATC [ref] from 20-22
    String delLoc = "chr1";
    int delLocStart = 20;
    int delLocStop = 22;

    // - [ref] / ATC from 20-20
    String insLoc = "chr1";
    int insLocStart = 20;
    int insLocStop = 20;

    VariantContextBuilder basicBuilder, snpBuilder, insBuilder;

    @BeforeSuite
    public void before() {
        del = Allele.apply("A");
        delRef = Allele.apply("A", true);

        A = Allele.apply("A");
        C = Allele.apply("C");
        Aref = Allele.apply("A", true);
        T = Allele.apply("T");
        Tref = Allele.apply("T", true);

        ATC = Allele.apply("ATC");
        ATCref = Allele.apply("ATC", true);
    }

    @BeforeMethod
    public void beforeTest() {
        basicBuilder = new VariantContextBuilder("test", snpLoc,snpLocStart, snpLocStop, new AlleleContext(new Allele[]{Aref, T}));
        snpBuilder = new VariantContextBuilder("test", snpLoc,snpLocStart, snpLocStop, new AlleleContext(new Allele[]{Aref, T}));
        insBuilder = new VariantContextBuilder("test", insLoc, insLocStart, insLocStop, new AlleleContext(new Allele[]{delRef, ATC}));
    }

    @Test
    public void testDetermineTypes() {
        Allele ACref = Allele.apply("AC", true);
        Allele AC = Allele.apply("AC");
        Allele AT = Allele.apply("AT");
        Allele C = Allele.apply("C");
        Allele CAT = Allele.apply("CAT");
        Allele TAref = Allele.apply("TA", true);
        Allele TA = Allele.apply("TA");
        Allele TC = Allele.apply("TC");
        Allele symbolic = Allele.apply("<FOO>");

        // test REF
        Allele[] alleles = new Allele[]{Tref};
        VariantContext vc = snpBuilder.alleles(alleles).stop(snpLocStop).make();
        Assert.assertEquals(vc.alleleContext().getType(), VariantType.NO_VARIATION);

        // test SNPs
        alleles = new Allele[]{Tref, A};
        vc = snpBuilder.alleles(alleles).stop(snpLocStop).make();
        Assert.assertEquals(vc.alleleContext().getType(), VariantType.SNP);

        alleles = new Allele[]{Tref, A, C};
        vc = snpBuilder.alleles(alleles).stop(snpLocStop).make();
        Assert.assertEquals(vc.alleleContext().getType(), VariantType.SNP);

        // test MNPs
        alleles = new Allele[]{ACref, TA};
        vc = snpBuilder.alleles(alleles).stop(snpLocStop+1).make();
        Assert.assertEquals(vc.alleleContext().getType(), VariantType.MNP);

        alleles = new Allele[]{ATCref, CAT, Allele.apply("GGG")};
        vc = basicBuilder.alleles(alleles).stop(snpLocStop+2).make();
        Assert.assertEquals(vc.alleleContext().getType(), VariantType.MNP);

        // test INDELs
        alleles = new Allele[]{Aref, ATC};
        vc = basicBuilder.alleles(alleles).stop(snpLocStop).make();
        Assert.assertEquals(vc.alleleContext().getType(), VariantType.INDEL);

        alleles = new Allele[]{ATCref, A};
        vc = basicBuilder.alleles(alleles).stop(snpLocStop+2).make();
        Assert.assertEquals(vc.alleleContext().getType(), VariantType.INDEL);

        alleles = new Allele[]{Tref, TA, TC};
        vc = basicBuilder.alleles(alleles).stop(snpLocStop).make();
        Assert.assertEquals(vc.alleleContext().getType(), VariantType.INDEL);

        alleles = new Allele[]{ATCref, A, AC};
        vc = basicBuilder.alleles(alleles).stop(snpLocStop+2).make();
        Assert.assertEquals(vc.alleleContext().getType(), VariantType.INDEL);

        alleles = new Allele[]{ATCref, A, Allele.apply("ATCTC")};
        vc = basicBuilder.alleles(alleles).stop(snpLocStop+2).make();
        Assert.assertEquals(vc.alleleContext().getType(), VariantType.INDEL);

        // test MIXED
        alleles = new Allele[]{TAref, T, TC};
        vc = basicBuilder.alleles(alleles).stop(snpLocStop+1).make();
        Assert.assertEquals(vc.alleleContext().getType(), VariantType.MIXED);

        alleles = new Allele[]{TAref, T, AC};
        vc = basicBuilder.alleles(alleles).stop(snpLocStop+1).make();
        Assert.assertEquals(vc.alleleContext().getType(), VariantType.MIXED);

        alleles = new Allele[]{ACref, ATC, AT};
        vc = basicBuilder.alleles(alleles).stop(snpLocStop+1).make();
        Assert.assertEquals(vc.alleleContext().getType(), VariantType.MIXED);

        alleles = new Allele[]{Aref, T, symbolic};
        vc = basicBuilder.alleles(alleles).stop(snpLocStop).make();
        Assert.assertEquals(vc.alleleContext().getType(), VariantType.MIXED);

        // test SYMBOLIC
        alleles = new Allele[]{Tref, symbolic};
        vc = basicBuilder.alleles(alleles).stop(snpLocStop).make();
        Assert.assertEquals(vc.alleleContext().getType(), VariantType.SYMBOLIC);
    }

    @Test
    public void testMultipleSNPAlleleOrdering() {
        final Allele[] allelesNaturalOrder = new Allele[]{Aref, C, T};
        final Allele[]  allelesUnnaturalOrder = new Allele[]{Aref, T, C};
        VariantContext naturalVC = snpBuilder.alleles(allelesNaturalOrder).make();
        VariantContext unnaturalVC = snpBuilder.alleles(allelesUnnaturalOrder).make();
        Assert.assertTrue(Arrays.equals(naturalVC.alleleContext().getAlleles(), allelesNaturalOrder));
        Assert.assertTrue(Arrays.equals(unnaturalVC.alleleContext().getAlleles(), allelesUnnaturalOrder));
    }

    @Test
    public void testCreatingSNPVariantContext() {

        Allele[] alleles = new Allele[]{Aref, T};
        VariantContext vc = snpBuilder.alleles(alleles).make();

        Assert.assertEquals(vc.getChr(), snpLoc);
        Assert.assertEquals(vc.getStart(), snpLocStart);
        Assert.assertEquals(vc.getEnd(), snpLocStop);
        Assert.assertEquals(vc.alleleContext().getType(), VariantType.SNP);
        Assert.assertTrue(vc.alleleContext().isSNP());
        Assert.assertFalse(vc.alleleContext().isIndel());
        Assert.assertFalse(vc.alleleContext().isSimpleInsertion());
        Assert.assertFalse(vc.alleleContext().isSimpleDeletion());
        Assert.assertFalse(vc.alleleContext().isSimpleIndel());
        Assert.assertFalse(vc.alleleContext().isMixed());
        Assert.assertTrue(vc.alleleContext().isBiallelic());
        Assert.assertEquals(vc.alleleContext().getNAlleles(), 2);

        Assert.assertEquals(vc.alleleContext().getReference(), Aref);
        Assert.assertEquals(vc.alleleContext().getAlleles().length, 2);
        Assert.assertEquals(vc.alleleContext().getAlternateAlleles().length, 1);
        Assert.assertEquals(vc.alleleContext().getAlternateAllele(0), T);

        Assert.assertFalse(vc.hasGenotypes());

        Assert.assertEquals(vc.genotypeContext().getSampleNames().size(), 0);
    }

    @Test
    public void testCreatingRefVariantContext() {
        Allele[] alleles = new Allele[]{Aref};
        VariantContext vc = snpBuilder.alleles(alleles).make();

        Assert.assertEquals(vc.getChr(), snpLoc);
        Assert.assertEquals(vc.getStart(), snpLocStart);
        Assert.assertEquals(vc.getEnd(), snpLocStop);
        Assert.assertEquals(VariantType.NO_VARIATION, vc.alleleContext().getType());
        Assert.assertFalse(vc.alleleContext().isSNP());
        Assert.assertFalse(vc.alleleContext().isIndel());
        Assert.assertFalse(vc.alleleContext().isSimpleInsertion());
        Assert.assertFalse(vc.alleleContext().isSimpleDeletion());
        Assert.assertFalse(vc.alleleContext().isSimpleIndel());
        Assert.assertFalse(vc.alleleContext().isMixed());
        Assert.assertFalse(vc.alleleContext().isBiallelic());
        Assert.assertEquals(vc.alleleContext().getNAlleles(), 1);

        Assert.assertEquals(vc.alleleContext().getReference(), Aref);
        Assert.assertEquals(vc.alleleContext().getAlleles().length, 1);
        Assert.assertEquals(vc.alleleContext().getAlternateAlleles().length, 0);
        //Assert.assertEquals(vc.getAlternateAllele(0), T);

        Assert.assertFalse(vc.hasGenotypes());
        Assert.assertEquals(vc.genotypeContext().getSampleNames().size(), 0);
    }

    @Test
    public void testCreatingDeletionVariantContext() {
        Allele[] alleles = new Allele[]{ATCref, del};
        VariantContext vc = new VariantContextBuilder("test", delLoc, delLocStart, delLocStop, new AlleleContext(alleles)).make();

        Assert.assertEquals(vc.getChr(), delLoc);
        Assert.assertEquals(vc.getStart(), delLocStart);
        Assert.assertEquals(vc.getEnd(), delLocStop);
        Assert.assertEquals(vc.alleleContext().getType(), VariantType.INDEL);
        Assert.assertFalse(vc.alleleContext().isSNP());
        Assert.assertTrue(vc.alleleContext().isIndel());
        Assert.assertFalse(vc.alleleContext().isSimpleInsertion());
        Assert.assertTrue(vc.alleleContext().isSimpleDeletion());
        Assert.assertTrue(vc.alleleContext().isSimpleIndel());
        Assert.assertFalse(vc.alleleContext().isMixed());
        Assert.assertTrue(vc.alleleContext().isBiallelic());
        Assert.assertEquals(vc.alleleContext().getNAlleles(), 2);

        Assert.assertEquals(vc.alleleContext().getReference(), ATCref);
        Assert.assertEquals(vc.alleleContext().getAlleles().length, 2);
        Assert.assertEquals(vc.alleleContext().getAlternateAlleles().length, 1);
        Assert.assertEquals(vc.alleleContext().getAlternateAllele(0), del);

        Assert.assertFalse(vc.hasGenotypes());

        Assert.assertEquals(vc.genotypeContext().getSampleNames().size(), 0);
    }

    @Test
    public void testCreatingComplexSubstitutionVariantContext() {
        Allele[] alleles = new Allele[]{Tref, ATC};
        VariantContext vc = new VariantContextBuilder("test", insLoc, insLocStart, insLocStop, new AlleleContext(alleles)).make();

        Assert.assertEquals(vc.getChr(), insLoc);
        Assert.assertEquals(vc.getStart(), insLocStart);
        Assert.assertEquals(vc.getEnd(), insLocStop);
        Assert.assertEquals(vc.alleleContext().getType(), VariantType.INDEL);
        Assert.assertFalse(vc.alleleContext().isSNP());
        Assert.assertTrue(vc.alleleContext().isIndel());
        Assert.assertFalse(vc.alleleContext().isSimpleInsertion());
        Assert.assertFalse(vc.alleleContext().isSimpleDeletion());
        Assert.assertFalse(vc.alleleContext().isSimpleIndel());
        Assert.assertFalse(vc.alleleContext().isMixed());
        Assert.assertTrue(vc.alleleContext().isBiallelic());
        Assert.assertEquals(vc.alleleContext().getNAlleles(), 2);

        Assert.assertEquals(vc.alleleContext().getReference(), Tref);
        Assert.assertEquals(vc.alleleContext().getAlleles().length, 2);
        Assert.assertEquals(vc.alleleContext().getAlternateAlleles().length, 1);
        Assert.assertEquals(vc.alleleContext().getAlternateAllele(0), ATC);

        Assert.assertFalse(vc.hasGenotypes());

        Assert.assertEquals(vc.genotypeContext().getSampleNames().size(), 0);
    }

    @Test
    public void testMatchingAlleles() {
        Allele[] alleles = new Allele[]{ATCref, del};
        VariantContext vc = new VariantContextBuilder("test", delLoc, delLocStart, delLocStop, new AlleleContext(alleles)).make();
        VariantContext vc2 = new VariantContextBuilder("test2", delLoc, delLocStart+12, delLocStop+12, new AlleleContext(alleles)).make();

        Assert.assertTrue(vc.alleleContext().hasSameAllelesAs(vc2.alleleContext()));
        Assert.assertTrue(vc.alleleContext().hasSameAlternateAllelesAs(vc2.alleleContext()));
    }

    @Test
    public void testCreatingInsertionVariantContext() {
        Allele[] alleles = new Allele[]{delRef, ATC};
        VariantContext vc = insBuilder.alleles(alleles).make();

        Assert.assertEquals(vc.getChr(), insLoc);
        Assert.assertEquals(vc.getStart(), insLocStart);
        Assert.assertEquals(vc.getEnd(), insLocStop);
        Assert.assertEquals(vc.alleleContext().getType(), VariantType.INDEL);
        Assert.assertFalse(vc.alleleContext().isSNP());
        Assert.assertTrue(vc.alleleContext().isIndel());
        Assert.assertTrue(vc.alleleContext().isSimpleInsertion());
        Assert.assertFalse(vc.alleleContext().isSimpleDeletion());
        Assert.assertTrue(vc.alleleContext().isSimpleIndel());
        Assert.assertFalse(vc.alleleContext().isMixed());
        Assert.assertTrue(vc.alleleContext().isBiallelic());
        Assert.assertEquals(vc.alleleContext().getNAlleles(), 2);

        Assert.assertEquals(vc.alleleContext().getReference(), delRef);
        Assert.assertEquals(vc.alleleContext().getAlleles().length, 2);
        Assert.assertEquals(vc.alleleContext().getAlternateAlleles().length, 1);
        Assert.assertEquals(vc.alleleContext().getAlternateAllele(0), ATC);
        Assert.assertFalse(vc.hasGenotypes());

        Assert.assertEquals(vc.genotypeContext().getSampleNames().size(), 0);
    }

    @Test
    public void testCreatingPartiallyCalledGenotype() {
        Allele[] alleles = new Allele[]{Aref, C};
        Genotype g = GenotypeBuilder.apply("foo", new Allele[]{C, Allele.NO_CALL()});
        VariantContext vc = new VariantContextBuilder("test", snpLoc, snpLocStart, snpLocStop, new AlleleContext(alleles)).genotypes(scala.collection.JavaConversions.asScalaBuffer(Arrays.asList(g))).make();

        Assert.assertTrue(vc.alleleContext().isSNP());
        Assert.assertEquals(vc.alleleContext().getNAlleles(), 2);
        Assert.assertTrue(vc.hasGenotypes());
        Assert.assertFalse(vc.genotypeContext().isMonomorphicInSamples(Aref));
        Assert.assertTrue(vc.genotypeContext().isPolymorphicInSamples(Aref));
        Assert.assertEquals(vc.genotypeContext().apply("foo"), g);
        Assert.assertEquals(vc.genotypeContext().getCalledChrCount(), 1); // we only have 1 called chromosomes, we exclude the NO_CALL one isn't called
        Assert.assertEquals(vc.genotypeContext().getCalledChrCount(Aref), 0);
        Assert.assertEquals(vc.genotypeContext().getCalledChrCount(C), 1);
        Assert.assertFalse(vc.genotypeContext().apply("foo").isHet());
        Assert.assertFalse(vc.genotypeContext().apply("foo").isHom());
        Assert.assertFalse(vc.genotypeContext().apply("foo").isNoCall());
        Assert.assertFalse(vc.genotypeContext().apply("foo").isHom());
        Assert.assertTrue(vc.genotypeContext().apply("foo").isMixed());
        Assert.assertEquals(vc.genotypeContext().apply("foo").getType(), GenotypeType.MIXED);
    }

    @Test (expectedExceptions = Exception.class)
    public void testBadConstructorArgs1() {
        new VariantContextBuilder("test", insLoc, insLocStart, insLocStop, new AlleleContext(new Allele[]{delRef, ATCref})).make();
    }

    @Test (expectedExceptions = Exception.class)
    public void testBadConstructorArgs2() {
        new VariantContextBuilder("test", insLoc, insLocStart, insLocStop, new AlleleContext(new Allele[]{delRef, del})).make();
    }

    @Test (expectedExceptions = Exception.class)
    public void testBadConstructorArgs3() {
        new VariantContextBuilder("test", insLoc, insLocStart, insLocStop, new AlleleContext(new Allele[]{del})).make();
    }

    @Test (expectedExceptions = Throwable.class)
    public void testBadConstructorArgs4() {
        new VariantContextBuilder("test", insLoc, insLocStart, insLocStop, new AlleleContext(new Allele[]{})).make();
    }

    @Test (expectedExceptions = Exception.class)
    public void testBadConstructorArgsDuplicateAlleles1() {
        new VariantContextBuilder("test", insLoc, insLocStart, insLocStop, new AlleleContext(new Allele[]{Aref, T, T})).make();
    }

    @Test (expectedExceptions = Exception.class)
    public void testBadConstructorArgsDuplicateAlleles2() {
        new VariantContextBuilder("test", insLoc, insLocStart, insLocStop, new AlleleContext(new Allele[]{Aref, A})).make();
    }

    @Test (expectedExceptions = Throwable.class)
    public void testBadLoc1() {
        Allele[] alleles = new Allele[]{Aref, T, del};
        new VariantContextBuilder("test", delLoc, delLocStart, delLocStop, new AlleleContext(alleles)).make();
    }

    @Test (expectedExceptions = Throwable.class)
    public void testBadID1() {
        new VariantContextBuilder("test", delLoc, delLocStart, delLocStop, new AlleleContext(new Allele[]{Aref, T})).id(null).make();
    }

    @Test (expectedExceptions = Exception.class)
    public void testBadID2() {
        new VariantContextBuilder("test", delLoc, delLocStart, delLocStop, new AlleleContext(new Allele[]{Aref, T})).id("").make();
    }

    @Test (expectedExceptions = Throwable.class)
    public void testBadPError() {
        new VariantContextBuilder("test", insLoc, insLocStart, insLocStop, new AlleleContext(new Allele[]{delRef, ATCref})).log10PError(0.5).make();
    }

    @Test
    public void testAccessingSimpleSNPGenotypes() {
        Allele[] alleles = new Allele[]{Aref, T};

        Genotype g1 = GenotypeBuilder.apply("AA", new Allele[]{Aref, Aref});
        Genotype g2 = GenotypeBuilder.apply("AT", new Allele[]{Aref, T});
        Genotype g3 = GenotypeBuilder.apply("TT", new Allele[]{T, T});

        VariantContext vc = new VariantContextBuilder("test", snpLoc, snpLocStart, snpLocStop, new AlleleContext(alleles))
                .genotypes(scala.collection.JavaConversions.asScalaBuffer(Arrays.asList(g1, g2, g3))).make();

        Assert.assertTrue(vc.hasGenotypes());
        Assert.assertFalse(vc.genotypeContext().isMonomorphicInSamples(Aref));
        Assert.assertTrue(vc.genotypeContext().isPolymorphicInSamples(Aref));
        Assert.assertEquals(vc.genotypeContext().getSampleNames().size(), 3);

        Assert.assertEquals(vc.genotypeContext().size(), 3);
        Assert.assertEquals(vc.genotypeContext().apply("AA"), g1);
        Assert.assertEquals(vc.genotypeContext().apply("AT"), g2);
        Assert.assertEquals(vc.genotypeContext().apply("TT"), g3);


        Assert.assertTrue(vc.genotypeContext().containsSample("AA"));
        Assert.assertTrue(vc.genotypeContext().containsSample("AT"));
        Assert.assertTrue(vc.genotypeContext().containsSample("TT"));
        Assert.assertFalse(vc.genotypeContext().containsSample("foo"));
        Assert.assertFalse(vc.genotypeContext().containsSample("TTT"));
        Assert.assertFalse(vc.genotypeContext().containsSample("at"));
        Assert.assertFalse(vc.genotypeContext().containsSample("tt"));

        Assert.assertEquals(vc.genotypeContext().getCalledChrCount(), 6);
        Assert.assertEquals(vc.genotypeContext().getCalledChrCount(Aref), 3);
        Assert.assertEquals(vc.genotypeContext().getCalledChrCount(T), 3);
    }

    @Test
    public void testAccessingCompleteGenotypes() {
        Allele[] alleles = new Allele[]{Aref, T, ATC};

        Genotype g1 = GenotypeBuilder.apply("AA", new Allele[]{Aref, Aref});
        Genotype g2 = GenotypeBuilder.apply("AT", new Allele[]{Aref, T});
        Genotype g3 = GenotypeBuilder.apply("TT", new Allele[]{T, T});
        Genotype g4 = GenotypeBuilder.apply("Td", new Allele[]{T, ATC});
        Genotype g5 = GenotypeBuilder.apply("dd", new Allele[]{ATC, ATC});
        Genotype g6 = GenotypeBuilder.apply("..", new Allele[]{Allele.NO_CALL(), Allele.NO_CALL()});

        VariantContext vc = new VariantContextBuilder("test", snpLoc, snpLocStart, snpLocStop, new AlleleContext(alleles))
                .genotypes(scala.collection.JavaConversions.asScalaBuffer(Arrays.asList(g1, g2, g3, g4, g5, g6))).make();

        Assert.assertTrue(vc.hasGenotypes());
        Assert.assertFalse(vc.genotypeContext().isMonomorphicInSamples(Aref));
        Assert.assertTrue(vc.genotypeContext().isPolymorphicInSamples(Aref));
        Assert.assertEquals(vc.genotypeContext().size(), 6);

        Assert.assertEquals(3, vc.genotypeContext().getGenotypesByName(scala.collection.JavaConversions.asScalaBuffer(Arrays.asList("AA", "Td", "dd"))).size());

        Assert.assertEquals(10, vc.genotypeContext().getCalledChrCount());
        Assert.assertEquals(3, vc.genotypeContext().getCalledChrCount(Aref));
        Assert.assertEquals(4, vc.genotypeContext().getCalledChrCount(T));
        Assert.assertEquals(3, vc.genotypeContext().getCalledChrCount(ATC));
        Assert.assertEquals(2, vc.genotypeContext().getCalledChrCount(Allele.NO_CALL()));
    }

    @Test
    public void testAccessingRefGenotypes() {
        Allele[] alleles1 = new Allele[]{Aref, T};
        Allele[] alleles2 = new Allele[]{Aref};
        Allele[] alleles3 = new Allele[]{Aref, T};
        for ( Allele[] alleles : Arrays.asList(alleles1, alleles2, alleles3)) {
            Genotype g1 = GenotypeBuilder.apply("AA1", new Allele[]{Aref, Aref});
            Genotype g2 = GenotypeBuilder.apply("AA2", new Allele[]{Aref, Aref});
            Genotype g3 = GenotypeBuilder.apply("..", new Allele[]{Allele.NO_CALL(), Allele.NO_CALL()});
            VariantContext vc = new VariantContextBuilder("test", snpLoc, snpLocStart, snpLocStop, new AlleleContext(alleles))
                    .genotypes(scala.collection.JavaConversions.asScalaBuffer(Arrays.asList(g1, g2, g3))).make();

            Assert.assertTrue(vc.hasGenotypes());
            Assert.assertTrue(vc.genotypeContext().isMonomorphicInSamples(Aref));
            Assert.assertFalse(vc.genotypeContext().isPolymorphicInSamples(Aref));
            Assert.assertEquals(vc.genotypeContext().size(), 3);

            Assert.assertEquals(4, vc.genotypeContext().getCalledChrCount());
            Assert.assertEquals(4, vc.genotypeContext().getCalledChrCount(Aref));
            Assert.assertEquals(0, vc.genotypeContext().getCalledChrCount(T));
            Assert.assertEquals(2, vc.genotypeContext().getCalledChrCount(Allele.NO_CALL()));
        }
    }

    @Test
    public void testFilters() {
        Allele[] alleles = new Allele[]{Aref, T};
        Genotype g1 = GenotypeBuilder.apply("AA", new Allele[]{Aref, Aref});
        Genotype g2 = GenotypeBuilder.apply("AT", new Allele[]{Aref, T});

        VariantContext vc = new VariantContextBuilder("test", snpLoc, snpLocStart, snpLocStop, new AlleleContext(alleles)).genotypes(scala.collection.JavaConversions.asScalaBuffer(Arrays.asList(g1, g2))).make();

        Assert.assertTrue(vc.isNotFiltered());
        Assert.assertFalse(vc.isFiltered());
        Assert.assertEquals(0, vc.getFilters().size());
        Assert.assertFalse(vc.filtersWereApplied());
        Assert.assertNull(vc.getFiltersMaybeNull());

        vc = new VariantContextBuilder(vc).filters(scala.collection.JavaConversions.asScalaBuffer(Arrays.asList("BAD_SNP_BAD!"))).make();

        Assert.assertFalse(vc.isNotFiltered());
        Assert.assertTrue(vc.isFiltered());
        Assert.assertEquals(1, vc.getFilters().size());
        Assert.assertTrue(vc.filtersWereApplied());
        Assert.assertNotNull(vc.getFiltersMaybeNull());

        Set<String> filters = new HashSet<String>(Arrays.asList("BAD_SNP_BAD!", "REALLY_BAD_SNP", "CHRIST_THIS_IS_TERRIBLE"));
        vc = new VariantContextBuilder(vc).filters(scala.collection.JavaConversions.asScalaSet(filters)).make();

        Assert.assertFalse(vc.isNotFiltered());
        Assert.assertTrue(vc.isFiltered());
        Assert.assertEquals(3, vc.getFilters().size());
        Assert.assertTrue(vc.filtersWereApplied());
        Assert.assertNotNull(vc.getFiltersMaybeNull());
    }

    @Test
    public void testGetGenotypeCounts() {
        Allele[] alleles = new Allele[]{Aref, T};
        Genotype g1 = GenotypeBuilder.apply("AA", new Allele[]{Aref, Aref});
        Genotype g2 = GenotypeBuilder.apply("AT", new Allele[]{Aref, T});
        Genotype g3 = GenotypeBuilder.apply("TT", new Allele[]{T, T});
        Genotype g4 = GenotypeBuilder.apply("A.", new Allele[]{Aref, Allele.NO_CALL()});
        Genotype g5 = GenotypeBuilder.apply("..", new Allele[]{Allele.NO_CALL(), Allele.NO_CALL()});

        // we need to create a new VariantContext each time
        VariantContext vc = new VariantContextBuilder("foo", snpLoc, snpLocStart, snpLocStop, new AlleleContext(alleles)).genotypes(scala.collection.JavaConversions.asScalaBuffer(Arrays.asList(g1,g2,g3,g4,g5))).make();
        Assert.assertEquals(1, vc.genotypeContext().getHetCount());
        vc = new VariantContextBuilder("foo", snpLoc, snpLocStart, snpLocStop, new AlleleContext(alleles)).genotypes(scala.collection.JavaConversions.asScalaBuffer(Arrays.asList(g1,g2,g3,g4,g5))).make();
        Assert.assertEquals(1, vc.genotypeContext().getHomRefCount());
        vc = new VariantContextBuilder("foo", snpLoc, snpLocStart, snpLocStop, new AlleleContext(alleles)).genotypes(scala.collection.JavaConversions.asScalaBuffer(Arrays.asList(g1,g2,g3,g4,g5))).make();
        Assert.assertEquals(1, vc.genotypeContext().getHomVarCount());
        vc = new VariantContextBuilder("foo", snpLoc, snpLocStart, snpLocStop, new AlleleContext(alleles)).genotypes(scala.collection.JavaConversions.asScalaBuffer(Arrays.asList(g1,g2,g3,g4,g5))).make();
        Assert.assertEquals(1, vc.genotypeContext().getMixedCount());
        vc = new VariantContextBuilder("foo", snpLoc, snpLocStart, snpLocStop, new AlleleContext(alleles)).genotypes(scala.collection.JavaConversions.asScalaBuffer(Arrays.asList(g1,g2,g3,g4,g5))).make();
        Assert.assertEquals(1, vc.genotypeContext().getNoCallCount());
    }

    @Test
    public void testVCFfromGenotypes() {
        Allele[] alleles = new Allele[]{Aref, C, T};
        Genotype g1 = GenotypeBuilder.apply("AA", new Allele[]{Aref, Aref});
        Genotype g2 = GenotypeBuilder.apply("AT", new Allele[]{Aref, T});
        Genotype g3 = GenotypeBuilder.apply("TT", new Allele[]{T, T});
        Genotype g4 = GenotypeBuilder.apply("..", new Allele[]{Allele.NO_CALL(), Allele.NO_CALL()});
        Genotype g5 = GenotypeBuilder.apply("AC", new Allele[]{Aref, C});
        VariantContext vc = new VariantContextBuilder("genotypes", snpLoc, snpLocStart, snpLocStop, new AlleleContext(alleles)).genotypes(scala.collection.JavaConversions.asScalaBuffer(Arrays.asList(g1,g2,g3,g4,g5))).make();

        VariantContext vc12 = vc.subContextFromSamples(scala.collection.JavaConversions.asScalaSet(new HashSet<String>(Arrays.asList(g1.getSampleName(), g2.getSampleName()))), true);
        VariantContext vc1 = vc.subContextFromSamples(scala.collection.JavaConversions.asScalaSet(new HashSet<String>(Arrays.asList(g1.getSampleName()))), true);
        VariantContext vc23 = vc.subContextFromSamples(scala.collection.JavaConversions.asScalaSet(new HashSet<String>(Arrays.asList(g2.getSampleName(), g3.getSampleName()))), true);
        VariantContext vc4 = vc.subContextFromSamples(scala.collection.JavaConversions.asScalaSet(new HashSet<String>(Arrays.asList(g4.getSampleName()))), true);
        VariantContext vc14 = vc.subContextFromSamples(scala.collection.JavaConversions.asScalaSet(new HashSet<String>(Arrays.asList(g1.getSampleName(), g4.getSampleName()))), true);
        VariantContext vc125 = vc.subContextFromSamples(scala.collection.JavaConversions.asScalaSet(new HashSet<String>(Arrays.asList(g1.getSampleName(), g2.getSampleName(), g5.getSampleName()))), true);

        Assert.assertTrue(vc12.genotypeContext().isPolymorphicInSamples(Aref));
        Assert.assertTrue(vc23.genotypeContext().isPolymorphicInSamples(Aref));
        Assert.assertTrue(vc1.genotypeContext().isMonomorphicInSamples(Aref));
        Assert.assertTrue(vc4.genotypeContext().isMonomorphicInSamples(Aref));
        Assert.assertTrue(vc14.genotypeContext().isMonomorphicInSamples(Aref));
        Assert.assertTrue(vc125.genotypeContext().isPolymorphicInSamples(Aref));

        Assert.assertTrue(vc12.alleleContext().isSNP());
        Assert.assertTrue(vc12.alleleContext().isVariant());
        Assert.assertTrue(vc12.alleleContext().isBiallelic());

        Assert.assertFalse(vc1.alleleContext().isSNP());
        Assert.assertFalse(vc1.alleleContext().isVariant());
        Assert.assertFalse(vc1.alleleContext().isBiallelic());

        Assert.assertTrue(vc23.alleleContext().isSNP());
        Assert.assertTrue(vc23.alleleContext().isVariant());
        Assert.assertTrue(vc23.alleleContext().isBiallelic());

        Assert.assertFalse(vc4.alleleContext().isSNP());
        Assert.assertFalse(vc4.alleleContext().isVariant());
        Assert.assertFalse(vc4.alleleContext().isBiallelic());

        Assert.assertFalse(vc14.alleleContext().isSNP());
        Assert.assertFalse(vc14.alleleContext().isVariant());
        Assert.assertFalse(vc14.alleleContext().isBiallelic());

        Assert.assertTrue(vc125.alleleContext().isSNP());
        Assert.assertTrue(vc125.alleleContext().isVariant());
        Assert.assertFalse(vc125.alleleContext().isBiallelic());

        Assert.assertEquals(3, vc12.genotypeContext().getCalledChrCount(Aref));
        Assert.assertEquals(1, vc23.genotypeContext().getCalledChrCount(Aref));
        Assert.assertEquals(2, vc1.genotypeContext().getCalledChrCount(Aref));
        Assert.assertEquals(0, vc4.genotypeContext().getCalledChrCount(Aref));
        Assert.assertEquals(2, vc14.genotypeContext().getCalledChrCount(Aref));
        Assert.assertEquals(4, vc125.genotypeContext().getCalledChrCount(Aref));
    }

    public void testGetGenotypeMethods() {
        Genotype g1 = GenotypeBuilder.apply("AA", new Allele[]{Aref, Aref});
        Genotype g2 = GenotypeBuilder.apply("AT", new Allele[]{Aref, T});
        Genotype g3 = GenotypeBuilder.apply("TT", new Allele[]{T, T});
        GenotypesContext gc = GenotypesContext.create(scala.collection.JavaConversions.asScalaBuffer(Arrays.asList(g1, g2, g3)));
        VariantContext vc = new VariantContextBuilder("genotypes", snpLoc, snpLocStart, snpLocStop, new AlleleContext(new Allele[]{Aref, T})).genotypes(gc).make();

        Assert.assertEquals(vc.genotypeContext().apply("AA"), g1);
        Assert.assertEquals(vc.genotypeContext().apply("AT"), g2);
        Assert.assertEquals(vc.genotypeContext().apply("TT"), g3);
        Assert.assertEquals(vc.genotypeContext().apply("CC"), null);

        Assert.assertEquals(vc.genotypeContext(), gc);
        Assert.assertTrue(vc.genotypeContext().getGenotypesByName(scala.collection.JavaConversions.asScalaBuffer(Arrays.asList("AA", "AT"))).sameElements(scala.collection.JavaConversions.asScalaBuffer(Arrays.asList(g1, g2))));
        Assert.assertTrue(vc.genotypeContext().getGenotypesByName(scala.collection.JavaConversions.asScalaBuffer(Arrays.asList("AA", "TT"))).sameElements(scala.collection.JavaConversions.asScalaBuffer(Arrays.asList(g1, g3))));
        Assert.assertTrue(vc.genotypeContext().getGenotypesByName(scala.collection.JavaConversions.asScalaBuffer(Arrays.asList("AA", "AT", "TT"))).sameElements(scala.collection.JavaConversions.asScalaBuffer(Arrays.asList(g1, g2, g3))));
        Assert.assertTrue(vc.genotypeContext().getGenotypesByName(scala.collection.JavaConversions.asScalaBuffer(Arrays.asList("AA", "AT", "CC"))).sameElements(scala.collection.JavaConversions.asScalaBuffer(Arrays.asList(g1, g2))));


        Assert.assertEquals(vc.genotypeContext().apply(0), g1);
        Assert.assertEquals(vc.genotypeContext().apply(1), g2);
        Assert.assertEquals(vc.genotypeContext().apply(2), g3);
    }

    // --------------------------------------------------------------------------------
    //
    // Test allele merging
    //
    // --------------------------------------------------------------------------------

    private class GetAllelesTest {
        List<Allele> alleles;
        String name;

        private GetAllelesTest(String name, Allele... arg) {
            this.name = name;
            this.alleles = Arrays.asList(arg);
        }

        public String toString() {
            return String.format("%s input=%s", name, alleles);
        }
    }

    @DataProvider(name = "getAlleles")
    public Object[][] mergeAllelesData() {
        List<Object[]> tests = new ArrayList<Object[]>();

        tests.add(new Object[]{new GetAllelesTest("A*",   Aref)});
        tests.add(new Object[]{new GetAllelesTest("A*/C", Aref, C)});
        tests.add(new Object[]{new GetAllelesTest("A*/C/T", Aref, C, T)});
        tests.add(new Object[]{new GetAllelesTest("A*/T/C", Aref, T, C)});
        tests.add(new Object[]{new GetAllelesTest("A*/C/T/ATC", Aref, C, T, ATC)});
        tests.add(new Object[]{new GetAllelesTest("A*/T/C/ATC", Aref, T, C, ATC)});
        tests.add(new Object[]{new GetAllelesTest("A*/ATC/T/C", Aref, ATC, T, C)});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "getAlleles")
    public void testMergeAlleles(GetAllelesTest cfg) {
        final List<Allele> altAlleles = cfg.alleles.subList(1, cfg.alleles.size());
        final VariantContext vc = new VariantContextBuilder("test", snpLoc, snpLocStart, snpLocStop, new AlleleContext(cfg.alleles.toArray(new Allele[cfg.alleles.size()]))).make();

        Assert.assertEquals(vc.alleleContext().getAlleles(), cfg.alleles, "VC alleles not the same as input alleles");
        Assert.assertEquals(vc.alleleContext().getNAlleles(), cfg.alleles.size(), "VC getNAlleles not the same as input alleles size");
        Assert.assertEquals(vc.alleleContext().getAlternateAlleles(), altAlleles, "VC alt alleles not the same as input alt alleles");


        for ( int i = 0; i < cfg.alleles.size(); i++ ) {
            final Allele inputAllele = cfg.alleles.get(i);

            Assert.assertTrue(vc.alleleContext().hasAllele(inputAllele));
            if ( inputAllele.isReference() ) {
                final Allele nonRefVersion = Allele.apply(inputAllele.getBases(), false);
                Assert.assertTrue(vc.alleleContext().hasAllele(nonRefVersion, true));
                Assert.assertFalse(vc.alleleContext().hasAllele(nonRefVersion, false));
            }

            Assert.assertEquals(inputAllele, vc.alleleContext().getAllele(inputAllele.getBaseString()));
            Assert.assertEquals(inputAllele, vc.alleleContext().getAllele(inputAllele.getBases()));

            if ( i > 0 ) { // it's an alt allele
                Assert.assertEquals(inputAllele, vc.alleleContext().getAlternateAllele(i - 1));
            }
        }

        final Allele missingAllele = Allele.apply("AACCGGTT"); // does not exist
        Assert.assertNull(vc.alleleContext().getAllele(missingAllele.getBases()));
        Assert.assertFalse(vc.alleleContext().hasAllele(missingAllele));
        Assert.assertFalse(vc.alleleContext().hasAllele(missingAllele, true));
    }

    private class SitesAndGenotypesVC {
        VariantContext vc, copy;
        String name;

        private SitesAndGenotypesVC(String name, VariantContext original) {
            this.name = name;
            this.vc = original;
            this.copy = new VariantContextBuilder(original).make();
        }

        public String toString() {
            return String.format("%s input=%s", name, vc);
        }
    }

    @DataProvider(name = "SitesAndGenotypesVC")
    public Object[][] MakeSitesAndGenotypesVCs() {
        Genotype g1 = GenotypeBuilder.apply("AA", new Allele[]{Aref, Aref});
        Genotype g2 = GenotypeBuilder.apply("AT", new Allele[]{Aref, T});
        Genotype g3 = GenotypeBuilder.apply("TT", new Allele[]{T, T});

        VariantContext sites = new VariantContextBuilder("sites", snpLoc, snpLocStart, snpLocStop, new AlleleContext( new Allele[]{Aref, T})).make();
        VariantContext genotypes = new VariantContextBuilder(sites).source("genotypes").genotypes(scala.collection.JavaConversions.asScalaBuffer(Arrays.asList(g1, g2, g3))).make();

        List<Object[]> tests = new ArrayList<Object[]>();

        tests.add(new Object[]{new SitesAndGenotypesVC("sites", sites)});
        tests.add(new Object[]{new SitesAndGenotypesVC("genotypes", genotypes)});

        return tests.toArray(new Object[][]{});
    }

    // --------------------------------------------------------------------------------
    //
    // Test modifying routines
    //
    // --------------------------------------------------------------------------------
    @Test(dataProvider = "SitesAndGenotypesVC")
    public void runModifyVCTests(SitesAndGenotypesVC cfg) {
        VariantContext modified = new VariantContextBuilder(cfg.vc).loc("chr2", 123, 123).make();
        Assert.assertEquals(modified.getChr(), "chr2");
        Assert.assertEquals(modified.getStart(), 123);
        Assert.assertEquals(modified.getEnd(), 123);

        modified = new VariantContextBuilder(cfg.vc).id("newID").make();
        Assert.assertEquals(modified.getID(), "newID");

        Set<String> newFilters = Collections.singleton("newFilter");
        modified = new VariantContextBuilder(cfg.vc).filters(scala.collection.JavaConversions.asScalaSet(newFilters)).make();
        Assert.assertEquals(modified.getFilters(), newFilters);

        // test the behavior when the builder's attribute object is null
        modified = new VariantContextBuilder(modified).attributes(null).make();
        Assert.assertTrue(modified.getAttributes().isEmpty());
        modified = new VariantContextBuilder(modified).attributes(null).rmAttribute("AC").make();
        Assert.assertTrue(modified.getAttributes().isEmpty());
        modified = new VariantContextBuilder(modified).attributes(null).attribute("AC", 1).make();
        Assert.assertEquals(modified.getAttribute("AC"), 1);

        // test the behavior when the builder's attribute object is not initialized
        modified = new VariantContextBuilder(modified.getSource(), modified.getChr(), modified.getStart(), modified.getEnd(), new AlleleContext(modified.alleleContext().getAlleles())).attribute("AC", 1).make();

        // test normal attribute modification
        modified = new VariantContextBuilder(cfg.vc).attribute("AC", 1).make();
        Assert.assertEquals(modified.getAttribute("AC"), 1);
        modified = new VariantContextBuilder(modified).attribute("AC", 2).make();
        Assert.assertEquals(modified.getAttribute("AC"), 2);

        Genotype g1 = GenotypeBuilder.apply("AA2", new Allele[]{Aref, Aref});
        Genotype g2 = GenotypeBuilder.apply("AT2", new Allele[]{Aref, T});
        Genotype g3 = GenotypeBuilder.apply("TT2", new Allele[]{T, T});
        GenotypesContext gc = GenotypesContext.create(scala.collection.JavaConversions.asScalaBuffer(Arrays.asList(g1, g2, g3)));
        modified = new VariantContextBuilder(cfg.vc).genotypes(gc).make();
        Assert.assertEquals(modified.genotypeContext(), gc);
        modified = new VariantContextBuilder(cfg.vc).noGenotypes().make();
        Assert.assertTrue(modified.genotypeContext().isEmpty());

        // test that original hasn't changed
        Assert.assertEquals(cfg.vc.getChr(), cfg.copy.getChr());
        Assert.assertEquals(cfg.vc.getStart(), cfg.copy.getStart());
        Assert.assertEquals(cfg.vc.getEnd(), cfg.copy.getEnd());
        Assert.assertEquals(cfg.vc.alleleContext().getAlleles(), cfg.copy.alleleContext().getAlleles());
        Assert.assertEquals(cfg.vc.getAttributes(), cfg.copy.getAttributes());
        Assert.assertEquals(cfg.vc.getID(), cfg.copy.getID());
        Assert.assertEquals(cfg.vc.genotypeContext(), cfg.copy.genotypeContext());
        Assert.assertEquals(cfg.vc.getLog10PError(), cfg.copy.getLog10PError());
        Assert.assertEquals(cfg.vc.getFilters(), cfg.copy.getFilters());
    }

    // --------------------------------------------------------------------------------
    //
    // Test subcontext
    //
    // --------------------------------------------------------------------------------
    private class SubContextTest {
        Set<String> samples;
        boolean updateAlleles;

        private SubContextTest(Collection<String> samples, boolean updateAlleles) {
            this.samples = new HashSet<String>(samples);
            this.updateAlleles = updateAlleles;
        }

        public String toString() {
            return String.format("%s samples=%s updateAlleles=%b", "SubContextTest", samples, updateAlleles);
        }
    }

    @DataProvider(name = "SubContextTest")
    public Object[][] MakeSubContextTest() {
        List<Object[]> tests = new ArrayList<Object[]>();

        for ( boolean updateAlleles : Arrays.asList(true, false)) {
            tests.add(new Object[]{new SubContextTest(Collections.<String>emptySet(), updateAlleles)});
            tests.add(new Object[]{new SubContextTest(Collections.singleton("MISSING"), updateAlleles)});
            tests.add(new Object[]{new SubContextTest(Collections.singleton("AA"), updateAlleles)});
            tests.add(new Object[]{new SubContextTest(Collections.singleton("AT"), updateAlleles)});
            tests.add(new Object[]{new SubContextTest(Collections.singleton("TT"), updateAlleles)});
            tests.add(new Object[]{new SubContextTest(Arrays.asList("AA", "AT"), updateAlleles)});
            tests.add(new Object[]{new SubContextTest(Arrays.asList("AA", "AT", "TT"), updateAlleles)});
            tests.add(new Object[]{new SubContextTest(Arrays.asList("AA", "AT", "MISSING"), updateAlleles)});
            tests.add(new Object[]{new SubContextTest(Arrays.asList("AA", "AT", "TT", "MISSING"), updateAlleles)});
            tests.add(new Object[]{new SubContextTest(Arrays.asList("AA", "AT", "AC"), updateAlleles)});
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "SubContextTest")
    public void runSubContextTest(SubContextTest cfg) {
        Genotype g1 = GenotypeBuilder.apply("AA", new Allele[]{Aref, Aref});
        Genotype g2 = GenotypeBuilder.apply("AT", new Allele[]{Aref, T});
        Genotype g3 = GenotypeBuilder.apply("TT", new Allele[]{T, T});
        Genotype g4 = GenotypeBuilder.apply("AC", new Allele[]{Aref, C});

        GenotypesContext gc = GenotypesContext.create(scala.collection.JavaConversions.asScalaBuffer(Arrays.asList(g1, g2, g3, g4)));
        VariantContext vc = new VariantContextBuilder("genotypes", snpLoc, snpLocStart, snpLocStop, new AlleleContext( new Allele[]{Aref, C, T})).genotypes(gc).make();
        VariantContext sub = vc.subContextFromSamples(scala.collection.JavaConversions.asScalaSet(cfg.samples), cfg.updateAlleles);

        // unchanged attributes should be the same
        Assert.assertEquals(sub.getChr(), vc.getChr());
        Assert.assertEquals(sub.getStart(), vc.getStart());
        Assert.assertEquals(sub.getEnd(), vc.getEnd());
        Assert.assertEquals(sub.getLog10PError(), vc.getLog10PError());
        Assert.assertEquals(sub.getFilters(), vc.getFilters());
        Assert.assertEquals(sub.getID(), vc.getID());
        Assert.assertEquals(sub.getAttributes(), vc.getAttributes());

        Set<Genotype> expectedGenotypes = new HashSet<Genotype>();
        if ( cfg.samples.contains(g1.getSampleName()) ) expectedGenotypes.add(g1);
        if ( cfg.samples.contains(g2.getSampleName()) ) expectedGenotypes.add(g2);
        if ( cfg.samples.contains(g3.getSampleName()) ) expectedGenotypes.add(g3);
        if ( cfg.samples.contains(g4.getSampleName()) ) expectedGenotypes.add(g4);
        GenotypesContext expectedGC = GenotypesContext.copy(scala.collection.JavaConversions.asScalaSet(expectedGenotypes));

        // these values depend on the results of sub
        if ( cfg.updateAlleles ) {
            // do the work to see what alleles should be here, and which not
            List<Allele> expectedAlleles = new ArrayList<Allele>();
            expectedAlleles.add(Aref);

            Set<Allele> genotypeAlleles = new HashSet<Allele>();
            for ( final Genotype g : expectedGC.getGenotypes() )
            {
                genotypeAlleles.addAll(Arrays.asList(g.getAlleles()));
            }
            genotypeAlleles.remove(Aref);

            // ensure original allele order
            for (Allele allele: vc.alleleContext().getAlleles())
                if (genotypeAlleles.contains(allele))
                    expectedAlleles.add(allele);

            Assert.assertEquals(sub.alleleContext().getAlleles(), expectedAlleles);
        } else {
            // not updating alleles -- should be the same
            Assert.assertEquals(sub.alleleContext().getAlleles(), vc.alleleContext().getAlleles());
        }

        // same sample names => success
        Assert.assertTrue(sub.genotypeContext().getSampleNames().equals(expectedGC.getSampleNames()));
    }

    // --------------------------------------------------------------------------------
    //
    // Test sample name functions
    //
    // --------------------------------------------------------------------------------
    private class SampleNamesTest {
        List<String> sampleNames;
        List<String> sampleNamesInOrder;

        private SampleNamesTest(List<String> sampleNames, List<String> sampleNamesInOrder) {
            this.sampleNamesInOrder = sampleNamesInOrder;
            this.sampleNames = sampleNames;
        }

        public String toString() {
            return String.format("%s samples=%s order=%s", "SampleNamesTest", sampleNames, sampleNamesInOrder);
        }
    }

    @DataProvider(name = "SampleNamesTest")
    public Object[][] MakeSampleNamesTest() {
        List<Object[]> tests = new ArrayList<Object[]>();

        tests.add(new Object[]{new SampleNamesTest(Arrays.asList("1"), Arrays.asList("1"))});
        tests.add(new Object[]{new SampleNamesTest(Arrays.asList("2", "1"), Arrays.asList("1", "2"))});
        tests.add(new Object[]{new SampleNamesTest(Arrays.asList("1", "2"), Arrays.asList("1", "2"))});
        tests.add(new Object[]{new SampleNamesTest(Arrays.asList("1", "2", "3"), Arrays.asList("1", "2", "3"))});
        tests.add(new Object[]{new SampleNamesTest(Arrays.asList("2", "1", "3"), Arrays.asList("1", "2", "3"))});
        tests.add(new Object[]{new SampleNamesTest(Arrays.asList("2", "3", "1"), Arrays.asList("1", "2", "3"))});
        tests.add(new Object[]{new SampleNamesTest(Arrays.asList("3", "1", "2"), Arrays.asList("1", "2", "3"))});
        tests.add(new Object[]{new SampleNamesTest(Arrays.asList("3", "2", "1"), Arrays.asList("1", "2", "3"))});
        tests.add(new Object[]{new SampleNamesTest(Arrays.asList("NA2", "NA1"), Arrays.asList("NA1", "NA2"))});

        return tests.toArray(new Object[][]{});
    }

    private final static void assertGenotypesAreInOrder(Iterable<Genotype> gIt, List<String> names) {
        int i = 0;
        for ( final Genotype g : gIt ) {
            Assert.assertEquals(g.getSampleName(), names.get(i), "Unexpected genotype ordering");
            i++;
        }
    }


    @Test(dataProvider = "SampleNamesTest")
    public void runSampleNamesTest(SampleNamesTest cfg) {
        GenotypesContext gc = GenotypesContext.create(cfg.sampleNames.size());
        for ( final String name : cfg.sampleNames ) {
            gc.$plus$eq(GenotypeBuilder.apply(name, new Allele[]{Aref, T}));
        }

        VariantContext vc = new VariantContextBuilder("genotypes", snpLoc, snpLocStart, snpLocStop, new AlleleContext(new Allele[]{Aref, T})).genotypes(gc).make();

        // same sample names => success
        Assert.assertTrue(vc.genotypeContext().getSampleNames().equals(new HashSet<String>(cfg.sampleNames)), "vc.genotypeContext().getSampleNames() = " + vc.genotypeContext().getSampleNames());
        Assert.assertEquals(vc.genotypeContext().getSampleNamesOrderedByName(), cfg.sampleNamesInOrder, "vc.getSampleNamesOrderedByName() = " + vc.genotypeContext().getSampleNamesOrderedByName());

        assertGenotypesAreInOrder(scala.collection.JavaConversions.asJavaIterable(vc.genotypeContext().iterateInSampleNameOrder()), cfg.sampleNamesInOrder);
        assertGenotypesAreInOrder(scala.collection.JavaConversions.asJavaIterable(vc.genotypeContext().iterateInSampleNameOrder(scala.collection.JavaConversions.asScalaIterable(cfg.sampleNames))), cfg.sampleNames);
    }

    @Test
    public void testGenotypeCounting() {
        Genotype noCall = GenotypeBuilder.apply("nocall", new Allele[]{Allele.NO_CALL()});
        Genotype mixed  = GenotypeBuilder.apply("mixed", new Allele[]{Aref, Allele.NO_CALL()});
        Genotype homRef = GenotypeBuilder.apply("homRef", new Allele[]{Aref, Aref});
        Genotype het    = GenotypeBuilder.apply("het", new Allele[]{Aref, T});
        Genotype homVar = GenotypeBuilder.apply("homVar", new Allele[]{T, T});

        List<Genotype> allGenotypes = Arrays.asList(noCall, mixed, homRef, het, homVar);
        final int nCycles = allGenotypes.size() * 10;

        for ( int i = 0; i < nCycles; i++ ) {
            int nNoCall = 0, nNoCallAlleles = 0, nA = 0, nT = 0, nMixed = 0, nHomRef = 0, nHet = 0, nHomVar = 0;
            int nSamples = 0;
            GenotypesContext gc = GenotypesContext.create();
            for ( int j = 0; j < i; j++ ) {
                nSamples++;
                Genotype g = allGenotypes.get(j % allGenotypes.size());
                final String name = String.format("%s_%d%d", g.getSampleName(), i, j);
                Genotype gt = GenotypeBuilder.apply(name, g.getAlleles());
                gc.$plus$eq(gt);
                switch ( g.getType() ) {
                    case NO_CALL: nNoCall++; nNoCallAlleles++; break;
                    case HOM_REF: nA += 2; nHomRef++; break;
                    case HET: nA++; nT++; nHet++; break;
                    case HOM_VAR: nT += 2; nHomVar++; break;
                    case MIXED: nA++; nNoCallAlleles++; nMixed++; break;
                    default: throw new RuntimeException("Unexpected genotype type " + g.getType());
                }

            }

            VariantContext vc = new VariantContextBuilder("genotypes", snpLoc, snpLocStart, snpLocStop, new AlleleContext(new Allele[]{Aref, T})).genotypes(gc).make();
            Assert.assertEquals(vc.genotypeContext().size(), nSamples);
            if ( nSamples > 0 ) {
                Assert.assertEquals(vc.genotypeContext().isPolymorphicInSamples(Aref), nT > 0);
                Assert.assertEquals(vc.genotypeContext().isMonomorphicInSamples(Aref), nT == 0);
            }
            Assert.assertEquals(vc.genotypeContext().getCalledChrCount(), nA + nT);

            Assert.assertEquals(vc.genotypeContext().getCalledChrCount(Allele.NO_CALL()), nNoCallAlleles);
            Assert.assertEquals(vc.genotypeContext().getCalledChrCount(Aref), nA);
            Assert.assertEquals(vc.genotypeContext().getCalledChrCount(T), nT);

            Assert.assertEquals(vc.genotypeContext().getNoCallCount(), nNoCall);
            Assert.assertEquals(vc.genotypeContext().getHomRefCount(), nHomRef);
            Assert.assertEquals(vc.genotypeContext().getHetCount(), nHet);
            Assert.assertEquals(vc.genotypeContext().getHomVarCount(), nHomVar);
            Assert.assertEquals(vc.genotypeContext().getMixedCount(), nMixed);
        }
    }
}

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

import org.testng.Assert
import org.testng.annotations.BeforeMethod
import org.testng.annotations.DataProvider
import org.testng.annotations.Test
import org.broadinstitute.variant.VariantBaseTest



class VariantContextUnitTest extends VariantBaseTest {

  val del = Allele("A")
  val delRef = Allele("A", true)

  val A = Allele("A")
  val C = Allele("C")
  val Aref = Allele("A", true)
  val T = Allele("T")
  val Tref = Allele("T", true)

  val ATC = Allele("ATC")
  val ATCref = Allele("ATC", true)


  // A [ref] / T at 10
  val snpLoc = "chr1"
  val snpLocStart = 10
  val snpLocStop = 10

  // - / ATC [ref] from 20-22
  val delLoc = "chr1"
  val delLocStart = 20
  val delLocStop = 22

  // - [ref] / ATC from 20-20
  val insLoc = "chr1"
  val insLocStart = 20
  val insLocStop = 20


  var basicBuilder : VariantContextBuilder = null
  var snpBuilder : VariantContextBuilder = null
  var insBuilder : VariantContextBuilder = null

  @BeforeMethod
  def beforeTest() {
    basicBuilder = new VariantContextBuilder("test", snpLoc,snpLocStart, snpLocStop, new AlleleContext(Array(Aref, T)))
    snpBuilder = new VariantContextBuilder("test", snpLoc,snpLocStart, snpLocStop, new AlleleContext(Array(Aref, T)))
    insBuilder = new VariantContextBuilder("test", insLoc, insLocStart, insLocStop, new AlleleContext(Array(delRef, ATC)))
  }




  @Test
  def testCreatingSNPVariantContext() {

    val alleles = Array(Aref, T)
    val vc = snpBuilder.alleles(alleles).make()

    Assert.assertEquals(vc.getChr(), snpLoc)
    Assert.assertEquals(vc.getStart(), snpLocStart)
    Assert.assertEquals(vc.getEnd(), snpLocStop)
    Assert.assertEquals(vc.alleleContext.getType(), VariantType.SNP)
    Assert.assertTrue(vc.alleleContext.isSNP())
    Assert.assertFalse(vc.alleleContext.isIndel())
    Assert.assertFalse(vc.alleleContext.isSimpleInsertion())
    Assert.assertFalse(vc.alleleContext.isSimpleDeletion())
    Assert.assertFalse(vc.alleleContext.isSimpleIndel())
    Assert.assertFalse(vc.alleleContext.isMixed())
    Assert.assertTrue(vc.alleleContext.isBiallelic())
    Assert.assertEquals(vc.alleleContext.getNAlleles(), 2)

    Assert.assertEquals(vc.alleleContext.getReference(), Aref)
    Assert.assertEquals(vc.alleleContext.getAlleles().size, 2)
    Assert.assertEquals(vc.alleleContext.getAlternateAlleles().size, 1)
    Assert.assertEquals(vc.alleleContext.getAlternateAllele(0), T)

    Assert.assertFalse(vc.hasGenotypes)

    Assert.assertEquals(vc.genotypeContext.getSampleNames().size, 0)
  }

  @Test
  def testCreatingRefVariantContext() {
    val alleles = Array(Aref)
    val vc = snpBuilder.alleles(alleles).make()

    Assert.assertEquals(vc.getChr(), snpLoc)
    Assert.assertEquals(vc.getStart(), snpLocStart)
    Assert.assertEquals(vc.getEnd(), snpLocStop)
    Assert.assertEquals(VariantType.NO_VARIATION, vc.alleleContext.getType())
    Assert.assertFalse(vc.alleleContext.isSNP())
    Assert.assertFalse(vc.alleleContext.isIndel())
    Assert.assertFalse(vc.alleleContext.isSimpleInsertion())
    Assert.assertFalse(vc.alleleContext.isSimpleDeletion())
    Assert.assertFalse(vc.alleleContext.isSimpleIndel())
    Assert.assertFalse(vc.alleleContext.isMixed())
    Assert.assertFalse(vc.alleleContext.isBiallelic())
    Assert.assertEquals(vc.alleleContext.getNAlleles(), 1)

    Assert.assertEquals(vc.alleleContext.getReference(), Aref)
    Assert.assertEquals(vc.alleleContext.getAlleles().size, 1)
    Assert.assertEquals(vc.alleleContext.getAlternateAlleles().size, 0)
    //Assert.assertEquals(vc.getAlternateAllele(0), T);

    Assert.assertFalse(vc.hasGenotypes)
    Assert.assertEquals(vc.genotypeContext.getSampleNames().size, 0)
  }

  @Test
  def testCreatingDeletionVariantContext() {
    val alleles = Array(ATCref, del)
    val vc = new VariantContextBuilder("test", delLoc, delLocStart, delLocStop, new AlleleContext(alleles)).make()

    Assert.assertEquals(vc.getChr(), delLoc)
    Assert.assertEquals(vc.getStart(), delLocStart)
    Assert.assertEquals(vc.getEnd(), delLocStop)
    Assert.assertEquals(vc.alleleContext.getType(), VariantType.INDEL)
    Assert.assertFalse(vc.alleleContext.isSNP())
    Assert.assertTrue(vc.alleleContext.isIndel())
    Assert.assertFalse(vc.alleleContext.isSimpleInsertion())
    Assert.assertTrue(vc.alleleContext.isSimpleDeletion())
    Assert.assertTrue(vc.alleleContext.isSimpleIndel())
    Assert.assertFalse(vc.alleleContext.isMixed())
    Assert.assertTrue(vc.alleleContext.isBiallelic())
    Assert.assertEquals(vc.alleleContext.getNAlleles(), 2)

    Assert.assertEquals(vc.alleleContext.getReference(), ATCref)
    Assert.assertEquals(vc.alleleContext.getAlleles().size, 2)
    Assert.assertEquals(vc.alleleContext.getAlternateAlleles().size, 1)
    Assert.assertEquals(vc.alleleContext.getAlternateAllele(0), del)

    Assert.assertFalse(vc.hasGenotypes)

    Assert.assertEquals(vc.genotypeContext.getSampleNames().size, 0)
  }

  @Test
  def testCreatingComplexSubstitutionVariantContext() {
    val alleles = Array(Tref, ATC)
    val vc = new VariantContextBuilder("test", insLoc, insLocStart, insLocStop, new AlleleContext(alleles)).make()

    Assert.assertEquals(vc.getChr(), insLoc)
    Assert.assertEquals(vc.getStart(), insLocStart)
    Assert.assertEquals(vc.getEnd(), insLocStop)
    Assert.assertEquals(vc.alleleContext.getType(), VariantType.INDEL)
    Assert.assertFalse(vc.alleleContext.isSNP())
    Assert.assertTrue(vc.alleleContext.isIndel())
    Assert.assertFalse(vc.alleleContext.isSimpleInsertion())
    Assert.assertFalse(vc.alleleContext.isSimpleDeletion())
    Assert.assertFalse(vc.alleleContext.isSimpleIndel())
    Assert.assertFalse(vc.alleleContext.isMixed())
    Assert.assertTrue(vc.alleleContext.isBiallelic())
    Assert.assertEquals(vc.alleleContext.getNAlleles(), 2)

    Assert.assertEquals(vc.alleleContext.getReference(), Tref)
    Assert.assertEquals(vc.alleleContext.getAlleles().size, 2)
    Assert.assertEquals(vc.alleleContext.getAlternateAlleles().size, 1)
    Assert.assertEquals(vc.alleleContext.getAlternateAllele(0), ATC)

    Assert.assertFalse(vc.hasGenotypes)

    Assert.assertEquals(vc.genotypeContext.getSampleNames().size, 0)
  }



  @Test
  def testCreatingInsertionVariantContext() {
    val alleles = Array(delRef, ATC)
    val vc = insBuilder.alleles(alleles).make()

    Assert.assertEquals(vc.getChr(), insLoc)
    Assert.assertEquals(vc.getStart(), insLocStart)
    Assert.assertEquals(vc.getEnd(), insLocStop)
    Assert.assertEquals(vc.alleleContext.getType(), VariantType.INDEL)
    Assert.assertFalse(vc.alleleContext.isSNP())
    Assert.assertTrue(vc.alleleContext.isIndel())
    Assert.assertTrue(vc.alleleContext.isSimpleInsertion())
    Assert.assertFalse(vc.alleleContext.isSimpleDeletion())
    Assert.assertTrue(vc.alleleContext.isSimpleIndel())
    Assert.assertFalse(vc.alleleContext.isMixed())
    Assert.assertTrue(vc.alleleContext.isBiallelic())
    Assert.assertEquals(vc.alleleContext.getNAlleles(), 2)

    Assert.assertEquals(vc.alleleContext.getReference(), delRef)
    Assert.assertEquals(vc.alleleContext.getAlleles().size, 2)
    Assert.assertEquals(vc.alleleContext.getAlternateAlleles().size, 1)
    Assert.assertEquals(vc.alleleContext.getAlternateAllele(0), ATC)
    Assert.assertFalse(vc.hasGenotypes)

    Assert.assertEquals(vc.genotypeContext.getSampleNames().size, 0)
  }

  @Test
  def testCreatingPartiallyCalledGenotype() {
    val alleles = Array(Aref, C)
    val g = GenotypeBuilder.apply("foo", Array(C, Allele.NO_CALL))
    val vc = new VariantContextBuilder("test", snpLoc, snpLocStart, snpLocStop, new AlleleContext(alleles)).genotypes(Array(g)).make()

    val reference = vc.alleleContext.getReference()

    Assert.assertTrue(vc.alleleContext.isSNP())
    Assert.assertEquals(vc.alleleContext.getNAlleles(), 2)
    Assert.assertTrue(vc.hasGenotypes)
    Assert.assertFalse(vc.genotypeContext.isMonomorphicInSamples(reference))
    Assert.assertTrue(vc.genotypeContext.isPolymorphicInSamples(reference))
    Assert.assertEquals(vc.genotypeContext("foo"), g)
    Assert.assertEquals(vc.genotypeContext.getCalledChrCount(), 1); // we only have 1 called chromosomes, we exclude the NO_CALL one isn't called
    Assert.assertEquals(vc.genotypeContext.getCalledChrCount(Aref), 0)
    Assert.assertEquals(vc.genotypeContext.getCalledChrCount(C), 1)
    Assert.assertFalse(vc.genotypeContext("foo").isHet())
    Assert.assertFalse(vc.genotypeContext("foo").isHom())
    Assert.assertFalse(vc.genotypeContext("foo").isNoCall())
    Assert.assertFalse(vc.genotypeContext("foo").isHom())
    Assert.assertTrue(vc.genotypeContext("foo").isMixed())
    Assert.assertEquals(vc.genotypeContext("foo").getType(), GenotypeType.MIXED)
  }

  @Test(expectedExceptions = Array(classOf[Exception]))
  def testBadConstructorArgs1() {
    new VariantContextBuilder("test", insLoc, insLocStart, insLocStop, new AlleleContext(Array(delRef, ATCref))).make()
  }

  @Test(expectedExceptions = Array(classOf[Exception]))
  def testBadConstructorArgs2() {
    new VariantContextBuilder("test", insLoc, insLocStart, insLocStop, new AlleleContext(Array(delRef, del))).make()
  }

  @Test(expectedExceptions = Array(classOf[Exception]))
  def testBadConstructorArgs3() {
    new VariantContextBuilder("test", insLoc, insLocStart, insLocStop, new AlleleContext(Array(del))).make()
  }

  @Test(expectedExceptions = Array(classOf[Exception]))
  def testBadConstructorArgs4() {
    new VariantContextBuilder("test", insLoc, insLocStart, insLocStop, new AlleleContext(Array[Allele]())).make()
  }

  @Test(expectedExceptions = Array(classOf[Exception]))
  def testBadConstructorArgsDuplicateAlleles1() {
    new VariantContextBuilder("test", insLoc, insLocStart, insLocStop, new AlleleContext(Array(Aref, T, T))).make()
  }

  @Test(expectedExceptions = Array(classOf[Exception]))
  def testBadConstructorArgsDuplicateAlleles2() {
    new VariantContextBuilder("test", insLoc, insLocStart, insLocStop, new AlleleContext(Array(Aref, A))).make()
  }

  @Test(expectedExceptions = Array(classOf[Exception]))
  def testBadLoc1() {
    val alleles = Array(Aref, T, del)
    new VariantContextBuilder("test", delLoc, delLocStart, delLocStop, new AlleleContext(alleles)).make()
  }

  @Test(expectedExceptions = Array(classOf[Exception]))
  def testBadID1() {
    new VariantContextBuilder("test", delLoc, delLocStart, delLocStop, new AlleleContext( Array(Aref, T))).id(null).make()
  }

  @Test(expectedExceptions = Array(classOf[Exception]))
  def testBadID2() {
    new VariantContextBuilder("test", delLoc, delLocStart, delLocStop, new AlleleContext(Array(Aref, T))).id("").make()
  }

  @Test(expectedExceptions = Array(classOf[Exception]))
  def testBadPError() {
    new VariantContextBuilder("test", insLoc, insLocStart, insLocStop, new AlleleContext(Array(delRef, ATCref))).log10PError(0.5).make()
  }

  @Test
  def testAccessingSimpleSNPGenotypes() {
    val alleles = Array(Aref, T)

    val g1 = GenotypeBuilder.apply("AA", Array(Aref, Aref))
    val g2 = GenotypeBuilder.apply("AT", Array(Aref, T))
    val g3 = GenotypeBuilder.apply("TT", Array(T, T))

    val vc = new VariantContextBuilder("test", snpLoc, snpLocStart, snpLocStop, new AlleleContext(alleles)).genotypes(Array(g1, g2, g3)).make()

    val reference = vc.alleleContext.getReference()

    Assert.assertTrue(vc.hasGenotypes)
    Assert.assertFalse(vc.genotypeContext.isMonomorphicInSamples(reference))
    Assert.assertTrue(vc.genotypeContext.isPolymorphicInSamples(reference))
    Assert.assertEquals(vc.genotypeContext.getSampleNames().size, 3)

    Assert.assertEquals(vc.genotypeContext.size, 3)
    Assert.assertEquals(vc.genotypeContext("AA"), g1)
    Assert.assertEquals(vc.genotypeContext("AT"), g2)
    Assert.assertEquals(vc.genotypeContext("TT"), g3)

    Assert.assertTrue(vc.genotypeContext.containsSample("AA"))
    Assert.assertTrue(vc.genotypeContext.containsSample("AT"))
    Assert.assertTrue(vc.genotypeContext.containsSample("TT"))
    Assert.assertFalse(vc.genotypeContext.containsSample("foo"))
    Assert.assertFalse(vc.genotypeContext.containsSample("TTT"))
    Assert.assertFalse(vc.genotypeContext.containsSample("at"))
    Assert.assertFalse(vc.genotypeContext.containsSample("tt"))

    Assert.assertEquals(vc.genotypeContext.getCalledChrCount(), 6)
    Assert.assertEquals(vc.genotypeContext.getCalledChrCount(Aref), 3)
    Assert.assertEquals(vc.genotypeContext.getCalledChrCount(T), 3)
  }

  @Test
  def testAccessingCompleteGenotypes() {
    val alleles = Array(Aref, T, ATC)

    val g1 = GenotypeBuilder.apply("AA", Array(Aref, Aref))
    val g2 = GenotypeBuilder.apply("AT", Array(Aref, T))
    val g3 = GenotypeBuilder.apply("TT", Array(T, T))
    val g4 = GenotypeBuilder.apply("Td", Array(T, ATC))
    val g5 = GenotypeBuilder.apply("dd", Array(ATC, ATC))
    val g6 = GenotypeBuilder.apply("..", Array(Allele.NO_CALL, Allele.NO_CALL))

    val vc = new VariantContextBuilder("test", snpLoc, snpLocStart, snpLocStop, new AlleleContext(alleles)).genotypes(Array(g1, g2, g3, g4, g5, g6)).make()

    Assert.assertTrue(vc.hasGenotypes)
    Assert.assertFalse(vc.genotypeContext.isMonomorphicInSamples(Aref))
    Assert.assertTrue(vc.genotypeContext.isPolymorphicInSamples(Aref))
    Assert.assertEquals(vc.genotypeContext.size, 6)

    Assert.assertEquals(3, vc.genotypeContext.subsetToSamples(Array("AA", "Td", "dd")).size)

    Assert.assertEquals(10, vc.genotypeContext.getCalledChrCount())
    Assert.assertEquals(3, vc.genotypeContext.getCalledChrCount(Aref))
    Assert.assertEquals(4, vc.genotypeContext.getCalledChrCount(T))
    Assert.assertEquals(3, vc.genotypeContext.getCalledChrCount(ATC))
    Assert.assertEquals(2, vc.genotypeContext.getCalledChrCount(Allele.NO_CALL))
  }

  @Test
  def testAccessingRefGenotypes() {
    val alleles1 = Array(Aref, T)
    val alleles2 = Array(Aref)
    val alleles3 = Array(Aref, T)
    for ( alleles <- Array(alleles1, alleles2, alleles3)) {
      val g1 = GenotypeBuilder.apply("AA1", Array(Aref, Aref))
      val g2 = GenotypeBuilder.apply("AA2", Array(Aref, Aref))
      val g3 = GenotypeBuilder.apply("..", Array(Allele.NO_CALL, Allele.NO_CALL))
      val vc = new VariantContextBuilder("test", snpLoc, snpLocStart, snpLocStop, new AlleleContext(alleles)).genotypes(Array(g1, g2, g3)).make()

      Assert.assertTrue(vc.hasGenotypes)
      Assert.assertTrue(vc.genotypeContext.isMonomorphicInSamples(Aref))
      Assert.assertFalse(vc.genotypeContext.isPolymorphicInSamples(Aref))
      Assert.assertEquals(vc.genotypeContext.size, 3)

      Assert.assertEquals(4, vc.genotypeContext.getCalledChrCount())
      Assert.assertEquals(4, vc.genotypeContext.getCalledChrCount(Aref))
      Assert.assertEquals(0, vc.genotypeContext.getCalledChrCount(T))
      Assert.assertEquals(2, vc.genotypeContext.getCalledChrCount(Allele.NO_CALL))
    }
  }

  @Test
  def testFilters() {
    val alleles = Array(Aref, T)
    val g1 = GenotypeBuilder.apply("AA", Array(Aref, Aref))
    val g2 = GenotypeBuilder.apply("AT", Array(Aref, T))

    var vc = new VariantContextBuilder("test", snpLoc, snpLocStart, snpLocStop, new AlleleContext(alleles)).genotypes(Array(g1, g2)).make()

    Assert.assertTrue(vc.isNotFiltered)
    Assert.assertFalse(vc.isFiltered)
    Assert.assertEquals(0, vc.getFilters.size)
    Assert.assertFalse(vc.filtersWereApplied)
    Assert.assertNull(vc.getFiltersMaybeNull)

    vc = new VariantContextBuilder(vc).filters("BAD_SNP_BAD!").make()

    Assert.assertFalse(vc.isNotFiltered)
    Assert.assertTrue(vc.isFiltered)
    Assert.assertEquals(1, vc.getFilters.size)
    Assert.assertTrue(vc.filtersWereApplied)
    Assert.assertNotNull(vc.getFiltersMaybeNull)

    val filters = new scala.collection.mutable.HashSet[String]() ++ Array("BAD_SNP_BAD!", "REALLY_BAD_SNP", "CHRIST_THIS_IS_TERRIBLE")
    vc = new VariantContextBuilder(vc).filters(filters).make()

    Assert.assertFalse(vc.isNotFiltered)
    Assert.assertTrue(vc.isFiltered)
    Assert.assertEquals(3, vc.getFilters.size)
    Assert.assertTrue(vc.filtersWereApplied)
    Assert.assertNotNull(vc.getFiltersMaybeNull)
  }

  @Test
  def testGetGenotypeCounts() {
    val alleles = Array(Aref, T)
    val g1 = GenotypeBuilder.apply("AA", Array(Aref, Aref))
    val g2 = GenotypeBuilder.apply("AT", Array(Aref, T))
    val g3 = GenotypeBuilder.apply("TT", Array(T, T))
    val g4 = GenotypeBuilder.apply("A.", Array(Aref, Allele.NO_CALL))
    val g5 = GenotypeBuilder.apply("..", Array(Allele.NO_CALL, Allele.NO_CALL))

    // we need to create a new VariantContext each time
    var vc = new VariantContextBuilder("foo", snpLoc, snpLocStart, snpLocStop, new AlleleContext(alleles)).genotypes(Array(g1,g2,g3,g4,g5)).make()
    Assert.assertEquals(1, vc.genotypeContext.getHetCount())
    vc = new VariantContextBuilder("foo", snpLoc, snpLocStart, snpLocStop, new AlleleContext(alleles)).genotypes(Array(g1,g2,g3,g4,g5)).make()
    Assert.assertEquals(1, vc.genotypeContext.getHomRefCount())
    vc = new VariantContextBuilder("foo", snpLoc, snpLocStart, snpLocStop, new AlleleContext(alleles)).genotypes(Array(g1,g2,g3,g4,g5)).make()
    Assert.assertEquals(1, vc.genotypeContext.getHomVarCount())
    vc = new VariantContextBuilder("foo", snpLoc, snpLocStart, snpLocStop, new AlleleContext(alleles)).genotypes(Array(g1,g2,g3,g4,g5)).make()
    Assert.assertEquals(1, vc.genotypeContext.getMixedCount())
    vc = new VariantContextBuilder("foo", snpLoc, snpLocStart, snpLocStop, new AlleleContext(alleles)).genotypes(Array(g1,g2,g3,g4,g5)).make()
    Assert.assertEquals(1, vc.genotypeContext.getNoCallCount())
  }

  @Test
  def testVCFfromGenotypes() {
    val alleles = Array(Aref, C, T)
    val g1 = GenotypeBuilder.apply("AA", Array(Aref, Aref))
    val g2 = GenotypeBuilder.apply("AT", Array(Aref, T))
    val g3 = GenotypeBuilder.apply("TT", Array(T, T))
    val g4 = GenotypeBuilder.apply("..", Array(Allele.NO_CALL, Allele.NO_CALL))
    val g5 = GenotypeBuilder.apply("AC", Array(Aref, C))
    val vc = new VariantContextBuilder("genotypes", snpLoc, snpLocStart, snpLocStop, new AlleleContext(alleles)).genotypes(Array(g1,g2,g3,g4,g5)).make()

    val vc12 = vc.subContextFromSamples(scala.collection.immutable.HashSet[String]() ++ Array(g1.getSampleName(), g2.getSampleName()), true)
    val vc1 = vc.subContextFromSamples(scala.collection.immutable.HashSet[String]() ++ Array(g1.getSampleName()), true)
    val vc23 = vc.subContextFromSamples(scala.collection.immutable.HashSet[String]() ++ Array(g2.getSampleName(), g3.getSampleName()), true)
    val vc4 = vc.subContextFromSamples(scala.collection.immutable.HashSet[String]() ++ Array(g4.getSampleName()), true)
    val vc14 = vc.subContextFromSamples(scala.collection.immutable.HashSet[String]() ++ Array(g1.getSampleName(), g4.getSampleName()), true)
    val vc125 = vc.subContextFromSamples(scala.collection.immutable.HashSet[String]() ++ Array(g1.getSampleName(), g2.getSampleName(), g5.getSampleName()), true)

    Assert.assertTrue(vc12.genotypeContext.isPolymorphicInSamples(Aref))
    Assert.assertTrue(vc23.genotypeContext.isPolymorphicInSamples(Aref))
    Assert.assertTrue(vc1.genotypeContext.isMonomorphicInSamples(Aref))
    Assert.assertTrue(vc4.genotypeContext.isMonomorphicInSamples(Aref))
    Assert.assertTrue(vc14.genotypeContext.isMonomorphicInSamples(Aref))
    Assert.assertTrue(vc125.genotypeContext.isPolymorphicInSamples(Aref))

    Assert.assertTrue(vc12.alleleContext.isSNP())
    Assert.assertTrue(vc12.alleleContext.isVariant())
    Assert.assertTrue(vc12.alleleContext.isBiallelic())

    Assert.assertFalse(vc1.alleleContext.isSNP())
    Assert.assertFalse(vc1.alleleContext.isVariant())
    Assert.assertFalse(vc1.alleleContext.isBiallelic())

    Assert.assertTrue(vc23.alleleContext.isSNP())
    Assert.assertTrue(vc23.alleleContext.isVariant())
    Assert.assertTrue(vc23.alleleContext.isBiallelic())

    Assert.assertFalse(vc4.alleleContext.isSNP())
    Assert.assertFalse(vc4.alleleContext.isVariant())
    Assert.assertFalse(vc4.alleleContext.isBiallelic())

    Assert.assertFalse(vc14.alleleContext.isSNP())
    Assert.assertFalse(vc14.alleleContext.isVariant())
    Assert.assertFalse(vc14.alleleContext.isBiallelic())

    Assert.assertTrue(vc125.alleleContext.isSNP())
    Assert.assertTrue(vc125.alleleContext.isVariant())
    Assert.assertFalse(vc125.alleleContext.isBiallelic())

    Assert.assertEquals(3, vc12.genotypeContext.getCalledChrCount(Aref))
    Assert.assertEquals(1, vc23.genotypeContext.getCalledChrCount(Aref))
    Assert.assertEquals(2, vc1.genotypeContext.getCalledChrCount(Aref))
    Assert.assertEquals(0, vc4.genotypeContext.getCalledChrCount(Aref))
    Assert.assertEquals(2, vc14.genotypeContext.getCalledChrCount(Aref))
    Assert.assertEquals(4, vc125.genotypeContext.getCalledChrCount(Aref))
  }

  def testGetGenotypeMethods() {
    val g1 = GenotypeBuilder.apply("AA", Array(Aref, Aref))
    val g2 = GenotypeBuilder.apply("AT", Array(Aref, T))
    val g3 = GenotypeBuilder.apply("TT", Array(T, T))
    val gc = GenotypesContext.create(g1, g2, g3)
    val vc = new VariantContextBuilder("genotypes", snpLoc, snpLocStart, snpLocStop, new AlleleContext(Array(Aref, T))).genotypeContext(gc).make()

    Assert.assertEquals(vc.genotypeContext("AA"), g1)
    Assert.assertEquals(vc.genotypeContext("AT"), g2)
    Assert.assertEquals(vc.genotypeContext("TT"), g3)
    Assert.assertEquals(vc.genotypeContext("CC"), null)

    Assert.assertEquals(vc.genotypeContext, gc)
    Assert.assertEquals(vc.genotypeContext.getGenotypesByName(Array("AA", "AT")), Array(g1, g2))
    Assert.assertEquals(vc.genotypeContext.getGenotypesByName(Array("AA", "TT")), Array(g1, g3))
    Assert.assertEquals(vc.genotypeContext.getGenotypesByName(Array("AA", "AT", "TT")), Array(g1, g2, g3))
    Assert.assertEquals(vc.genotypeContext.getGenotypesByName(Array("AA", "AT", "CC")), Array(g1, g2))

    Assert.assertEquals(vc.genotypeContext(0), g1)
    Assert.assertEquals(vc.genotypeContext(1), g2)
    Assert.assertEquals(vc.genotypeContext(2), g3)
  }

  // --------------------------------------------------------------------------------
  //
  // Test allele merging
  //
  // --------------------------------------------------------------------------------

  private class GetAllelesTest( val name : String, arg : Allele*) {
    val alleles = Array[Allele]() ++ arg

    override def toString() : String = {
      s"$name input=$alleles"
    }
  }

  @DataProvider(name = "getAlleles")
  def mergeAllelesData() : Array[Array[Any]] = {
    val tests = scala.collection.mutable.ArrayBuffer[Array[Any]]()

    tests += Array(new GetAllelesTest("A*",   Aref))
    tests += Array(new GetAllelesTest("A*/C", Aref, C))
    tests += Array(new GetAllelesTest("A*/C/T", Aref, C, T))
    tests += Array(new GetAllelesTest("A*/T/C", Aref, T, C))
    tests += Array(new GetAllelesTest("A*/C/T/ATC", Aref, C, T, ATC))
    tests += Array(new GetAllelesTest("A*/T/C/ATC", Aref, T, C, ATC))
    tests += Array(new GetAllelesTest("A*/ATC/T/C", Aref, ATC, T, C))

    tests.toArray
  }

  @Test(dataProvider = "getAlleles")
  def testMergeAlleles( cfg : GetAllelesTest) {
    val altAlleles = cfg.alleles.slice(1, cfg.alleles.size)
    val vc = new VariantContextBuilder("test", snpLoc, snpLocStart, snpLocStop, new AlleleContext(cfg.alleles)).make()

    Assert.assertEquals(vc.alleleContext.getAlleles(), cfg.alleles, "VC alleles not the same as input alleles")
    Assert.assertEquals(vc.alleleContext.getNAlleles(), cfg.alleles.size, "VC getNAlleles not the same as input alleles size")

    val altAllelesObserved = vc.alleleContext.getAlternateAlleles()
    Assert.assertTrue(altAllelesObserved.sameElements(altAlleles), "VC alt alleles not the same as input alt alleles")


    for ( i <- 0 until cfg.alleles.size ) {
      val inputAllele = cfg.alleles(i)

      Assert.assertTrue(vc.alleleContext.hasAllele(inputAllele))
      if ( inputAllele.isReference ) {
        val nonRefVersion = Allele(inputAllele.getBases(), false)
        Assert.assertTrue(vc.alleleContext.hasAllele(nonRefVersion, true))
        Assert.assertFalse(vc.alleleContext.hasAllele(nonRefVersion, false))
      }

      Assert.assertEquals(inputAllele, vc.alleleContext.getAllele(inputAllele.getBaseString()))
      Assert.assertEquals(inputAllele, vc.alleleContext.getAllele(inputAllele.getBases()))

      if ( i > 0 ) { // it's an alt allele
        Assert.assertEquals(inputAllele, vc.alleleContext.getAlternateAllele(i - 1))
      }
    }

    val missingAllele = Allele("AACCGGTT"); // does not exist
    Assert.assertNull(vc.alleleContext.getAllele(missingAllele.getBases()))
    Assert.assertFalse(vc.alleleContext.hasAllele(missingAllele))
    Assert.assertFalse(vc.alleleContext.hasAllele(missingAllele, true))
  }

  private class SitesAndGenotypesVC( val name : String,  original : VariantContext) {
    val vc = original
    val copy = new VariantContextBuilder(original).make()

    override def toString() : String =  {
      s"$name input=$vc"
    }
  }

  @DataProvider(name = "SitesAndGenotypesVC")
  def  MakeSitesAndGenotypesVCs() : Array[Array[Any]] = {
    val g1 = GenotypeBuilder.apply("AA", Array(Aref, Aref))
    val g2 = GenotypeBuilder.apply("AT", Array(Aref, T))
    val g3 = GenotypeBuilder.apply("TT", Array(T, T))

    val sites = new VariantContextBuilder("sites", snpLoc, snpLocStart, snpLocStop, new AlleleContext(Array(Aref, T))).make()
    val genotypes = new VariantContextBuilder(sites).source("genotypes").genotypes(Array(g1, g2, g3)).make()

    val tests = scala.collection.mutable.ArrayBuffer[Array[Any]]()

    tests += Array(new SitesAndGenotypesVC("sites", sites))
    tests += Array(new SitesAndGenotypesVC("genotypes", genotypes))

    tests.toArray
  }

  // --------------------------------------------------------------------------------
  //
  // Test modifying routines
  //
  // --------------------------------------------------------------------------------
  @Test(dataProvider = "SitesAndGenotypesVC")
  def runModifyVCTests( cfg : SitesAndGenotypesVC) {
    var modified = new VariantContextBuilder(cfg.vc).loc("chr2", 123, 123).make()
    Assert.assertEquals(modified.getChr(), "chr2")
    Assert.assertEquals(modified.getStart(), 123)
    Assert.assertEquals(modified.getEnd(), 123)

    modified = new VariantContextBuilder(cfg.vc).id("newID").make()
    Assert.assertEquals(modified.getID, "newID")

    val newFilters = scala.collection.immutable.HashSet() + "newFilter"
    modified = new VariantContextBuilder(cfg.vc).filters(newFilters).make()
    Assert.assertEquals(modified.getFilters, newFilters)

    // test the behavior when the builder's attribute object is null
    modified = new VariantContextBuilder(modified).attributes(null).make()
    Assert.assertTrue(modified.getAttributes.isEmpty)
    modified = new VariantContextBuilder(modified).attributes(null).rmAttribute("AC").make()
    Assert.assertTrue(modified.getAttributes.isEmpty)
    modified = new VariantContextBuilder(modified).attributes(null).attribute("AC", 1).make()
    Assert.assertEquals(modified.getAttribute("AC"), 1)

    // test the behavior when the builder's attribute object is not initialized
    modified = new VariantContextBuilder(modified.getSource, modified.getChr(), modified.getStart(), modified.getEnd(), modified.alleleContext).attribute("AC", 1).make();

    // test normal attribute modification
    modified = new VariantContextBuilder(cfg.vc).attribute("AC", 1).make()
    Assert.assertEquals(modified.getAttribute("AC"), 1)
    modified = new VariantContextBuilder(modified).attribute("AC", 2).make()
    Assert.assertEquals(modified.getAttribute("AC"), 2)

    val g1 = GenotypeBuilder.apply("AA2", Array(Aref, Aref))
    val g2 = GenotypeBuilder.apply("AT2", Array(Aref, T))
    val g3 = GenotypeBuilder.apply("TT2", Array(T, T))
    val gc = GenotypesContext.create(g1,g2,g3)
    modified = new VariantContextBuilder(cfg.vc).genotypes(gc).make()
    Assert.assertEquals(modified.genotypeContext, gc)
    modified = new VariantContextBuilder(cfg.vc).noGenotypes().make()
    Assert.assertTrue(modified.genotypeContext.isEmpty)

    // test that original hasn't changed
    Assert.assertEquals(cfg.vc.getChr(), cfg.copy.getChr())
    Assert.assertEquals(cfg.vc.getStart(), cfg.copy.getStart())
    Assert.assertEquals(cfg.vc.getEnd(), cfg.copy.getEnd())
    Assert.assertEquals(cfg.vc.alleleContext.getAlleles(), cfg.copy.alleleContext.getAlleles())
    Assert.assertEquals(cfg.vc.getAttributes, cfg.copy.getAttributes)
    Assert.assertEquals(cfg.vc.getID, cfg.copy.getID)
    Assert.assertEquals(cfg.vc.genotypeContext, cfg.copy.genotypeContext)
    Assert.assertEquals(cfg.vc.getLog10PError, cfg.copy.getLog10PError)
    Assert.assertEquals(cfg.vc.getFilters, cfg.copy.getFilters)
  }

  // --------------------------------------------------------------------------------
  //
  // Test subcontext
  //
  // --------------------------------------------------------------------------------
  private class SubContextTest(  _samples : Seq[String], val updateAlleles : Boolean) {
    val samples = _samples.toSet

    override def toString() : String = {
      s"SubContextTest samples=$samples updateAlleles=$updateAlleles"
    }
  }

  @DataProvider(name = "SubContextTest")
  def MakeSubContextTest() : Array[Array[Any]] = {
    val tests = scala.collection.mutable.ArrayBuffer[Array[Any]]()

    for ( updateAlleles  <- Array(true, false)) {
      tests += Array(new SubContextTest(Array[String](), updateAlleles))
      tests += Array(new SubContextTest(Array("MISSING"), updateAlleles))
      tests += Array(new SubContextTest(Array("AA"), updateAlleles))
      tests += Array(new SubContextTest(Array("AT"), updateAlleles))
      tests += Array(new SubContextTest(Array("TT"), updateAlleles))
      tests += Array(new SubContextTest(Array("AA", "AT"), updateAlleles))
      tests += Array(new SubContextTest(Array("AA", "AT", "TT"), updateAlleles))
      tests += Array(new SubContextTest(Array("AA", "AT", "MISSING"), updateAlleles))
      tests += Array(new SubContextTest(Array("AA", "AT", "TT", "MISSING"), updateAlleles))
      tests += Array(new SubContextTest(Array("AA", "AT", "AC"), updateAlleles))
    }

    tests.toArray
  }

  @Test(dataProvider = "SubContextTest")
  def runSubContextTest( cfg : SubContextTest) {
    val g1 = GenotypeBuilder.apply("AA", Array(Aref, Aref))
    val g2 = GenotypeBuilder.apply("AT", Array(Aref, T))
    val g3 = GenotypeBuilder.apply("TT", Array(T, T))
    val g4 = GenotypeBuilder.apply("AC", Array(Aref, C))

    val gc = GenotypesContext.create(g1, g2, g3, g4)
    val vc = new VariantContextBuilder("genotypes", snpLoc, snpLocStart, snpLocStop, new AlleleContext(Array(Aref, C, T))).genotypeContext(gc).make()
    val sub = vc.subContextFromSamples(cfg.samples, cfg.updateAlleles)

    // unchanged attributes should be the same
    Assert.assertEquals(sub.getChr(), vc.getChr())
    Assert.assertEquals(sub.getStart(), vc.getStart())
    Assert.assertEquals(sub.getEnd(), vc.getEnd())
    Assert.assertEquals(sub.getLog10PError, vc.getLog10PError)
    Assert.assertEquals(sub.getFilters, vc.getFilters)
    Assert.assertEquals(sub.getID, vc.getID)
    Assert.assertEquals(sub.getAttributes, vc.getAttributes)

    val expectedGenotypes = scala.collection.mutable.HashSet[Genotype]()
    if ( cfg.samples.contains(g1.getSampleName()) ) {expectedGenotypes += g1}
    if ( cfg.samples.contains(g2.getSampleName()) ) {expectedGenotypes += g2}
    if ( cfg.samples.contains(g3.getSampleName()) ) {expectedGenotypes += g3}
    if ( cfg.samples.contains(g4.getSampleName()) ) {expectedGenotypes += g4}
    val expectedGC = GenotypesContext.copy(expectedGenotypes)

    // these values depend on the results of sub
    if ( cfg.updateAlleles ) {
      // do the work to see what alleles should be here, and which not
      val expectedAlleles = scala.collection.mutable.ArrayBuffer[Allele]()
      expectedAlleles += Aref

      val genotypeAlleles = scala.collection.mutable.HashSet[Allele]()
      for (  g <- expectedGC )
        genotypeAlleles ++= g.getAlleles()
      genotypeAlleles -= Aref

      // ensure original allele order
      for ( allele <- vc.alleleContext.getAlleles())
        if (genotypeAlleles.contains(allele))
        {
          expectedAlleles += allele
        }

      Assert.assertTrue(sub.alleleContext.getAlleles().sameElements(expectedAlleles));
    } else {
      // not updating alleles -- should be the same
      Assert.assertEquals(sub.alleleContext.getAlleles(), vc.alleleContext.getAlleles());
    }

    // same sample names => success
    Assert.assertTrue(sub.genotypeContext.getSampleNames().equals(expectedGC.getSampleNames()));
  }

  // --------------------------------------------------------------------------------
  //
  // Test sample name functions
  //
  // --------------------------------------------------------------------------------
  private class SampleNamesTest( val sampleNames : Array[String],  val sampleNamesInOrder : Array[String]) {

    override def toString() : String = {
      s"SampleNamesTest samples=$sampleNames order=$sampleNamesInOrder"
    }
  }

  @DataProvider(name = "SampleNamesTest")
  def MakeSampleNamesTest() : Array[Array[Any]]  = {
    val tests = scala.collection.mutable.ArrayBuffer[Array[Any]]()

    tests += Array(new SampleNamesTest(Array("1"), Array("1")))
    tests += Array(new SampleNamesTest(Array("2", "1"), Array("1", "2")))
    tests += Array(new SampleNamesTest(Array("1", "2"), Array("1", "2")))
    tests += Array(new SampleNamesTest(Array("1", "2", "3"), Array("1", "2", "3")))
    tests += Array(new SampleNamesTest(Array("2", "1", "3"), Array("1", "2", "3")))
    tests += Array(new SampleNamesTest(Array("2", "3", "1"), Array("1", "2", "3")))
    tests += Array(new SampleNamesTest(Array("3", "1", "2"), Array("1", "2", "3")))
    tests += Array(new SampleNamesTest(Array("3", "2", "1"), Array("1", "2", "3")))
    tests += Array(new SampleNamesTest(Array("NA2", "NA1"), Array("NA1", "NA2")))

    tests.toArray
  }




  @Test(dataProvider = "SampleNamesTest")
  def runSampleNamesTest( cfg : SampleNamesTest) {
    val gc = GenotypesContext.create(cfg.sampleNames.size)
    for ( name <- cfg.sampleNames ) {
      gc += GenotypeBuilder.apply(name, Array(Aref, T))
    }

    val vc = new VariantContextBuilder("genotypes", snpLoc, snpLocStart, snpLocStop, new AlleleContext(Array(Aref, T))).genotypeContext(gc).make();

    // same sample names => success
    Assert.assertTrue(vc.genotypeContext.getSampleNames().equals( scala.collection.immutable.HashSet[String]() ++ cfg.sampleNames), "vc.getSampleNames() = " + vc.genotypeContext.getSampleNames());
    Assert.assertTrue(vc.genotypeContext.getSampleNamesOrderedByName().sameElements( cfg.sampleNamesInOrder), "vc.getSampleNamesOrderedByName() = " + vc.genotypeContext.getSampleNamesOrderedByName());

    VariantContextUnitTest.assertGenotypesAreInOrder(vc.genotypeContext.iterateInSampleNameOrder(), cfg.sampleNamesInOrder);
    VariantContextUnitTest.assertGenotypesAreInOrder(vc.genotypeContext.iterateInSampleNameOrder(cfg.sampleNames), cfg.sampleNames);
  }

  @Test
  def testGenotypeCounting() {
    val noCall = GenotypeBuilder.apply("nocall", Array(Allele.NO_CALL))
    val mixed  = GenotypeBuilder.apply("mixed", Array(Aref, Allele.NO_CALL))
    val homRef = GenotypeBuilder.apply("homRef", Array(Aref, Aref))
    val het    = GenotypeBuilder.apply("het", Array(Aref, T))
    val homVar = GenotypeBuilder.apply("homVar", Array(T, T))

    val allGenotypes = Array(noCall, mixed, homRef, het, homVar)
    val nCycles = allGenotypes.size * 10

    for ( i <- 0 until nCycles ) {
      var nNoCall = 0
      var nNoCallAlleles = 0
      var nA = 0
      var nT = 0
      var nMixed = 0
      var nHomRef = 0
      var nHet = 0
      var nHomVar = 0
      var nSamples = 0
      val gc = GenotypesContext.create();
      for (  j <- 0 until  i ) {
        nSamples += 1
        val g = allGenotypes(j % allGenotypes.size)
        val sampleName = g.getSampleName()
        val name = s"$sampleName$i$j"
        gc += GenotypeBuilder.apply(name, g.getAlleles())
        g.getType() match {
          case GenotypeType.NO_CALL =>  nNoCall += 1 ; nNoCallAlleles += 1
          case GenotypeType.HOM_REF =>  nA += 2; nHomRef += 1
          case GenotypeType.HET     =>  nA += 1; nT += 1 ; nHet += 1
          case GenotypeType.HOM_VAR =>  nT += 2; nHomVar += 1
          case GenotypeType.MIXED   =>  nA += 1; nNoCallAlleles += 1; nMixed += 1
          case _       => throw new RuntimeException("Unexpected genotype type " + g.getType());
        }

      }

      val vc = new VariantContextBuilder("genotypes", snpLoc, snpLocStart, snpLocStop, new AlleleContext(Array(Aref, T))).genotypeContext(gc).make();
      Assert.assertEquals(vc.genotypeContext.size, nSamples);
      if ( nSamples > 0 ) {
        Assert.assertEquals(vc.genotypeContext.isPolymorphicInSamples(Aref), nT > 0);
        Assert.assertEquals(vc.genotypeContext.isMonomorphicInSamples(Aref), nT == 0);
      }
      Assert.assertEquals(vc.genotypeContext.getCalledChrCount(), nA + nT);

      Assert.assertEquals(vc.genotypeContext.getCalledChrCount(Allele.NO_CALL), nNoCallAlleles);
      Assert.assertEquals(vc.genotypeContext.getCalledChrCount(Aref), nA);
      Assert.assertEquals(vc.genotypeContext.getCalledChrCount(T), nT);

      Assert.assertEquals(vc.genotypeContext.getNoCallCount(), nNoCall);
      Assert.assertEquals(vc.genotypeContext.getHomRefCount(), nHomRef);
      Assert.assertEquals(vc.genotypeContext.getHetCount(), nHet);
      Assert.assertEquals(vc.genotypeContext.getHomVarCount(), nHomVar);
      Assert.assertEquals(vc.genotypeContext.getMixedCount(), nMixed);
    }
  }
}

object VariantContextUnitTest
{
  private def assertGenotypesAreInOrder( gIt : Iterable[Genotype],  names : Array[String]) {
    var i = 0;
    for ( g <- gIt ) {
      Assert.assertEquals(g.getSampleName(), names(i), "Unexpected genotype ordering");
      i +=1
    }
  }

}

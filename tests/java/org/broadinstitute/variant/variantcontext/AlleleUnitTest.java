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
import org.testng.Assert;
import org.testng.annotations.BeforeSuite;
import org.testng.annotations.Test;

//    public Allele(byte[] bases, boolean isRef) {
//    public Allele(boolean isRef) {
//    public Allele(String bases, boolean isRef) {
//    public boolean isReference()        { return isRef; }
//    public boolean isNonReference()     { return ! isReference(); }
//    public byte[] getBases() { return bases; }
//    public boolean equals(Allele other) {
//    public int length() {

/**
 * Basic unit test for RecalData
 */
public class AlleleUnitTest extends VariantBaseTest {
    Allele ARef, A, T, ATIns, ATCIns, NoCall;
    
    @BeforeSuite
    public void before() {
        A = Allele.apply("A");
        ARef = Allele.apply("A", true);
        T = Allele.apply("T");

        ATIns = Allele.apply("AT");
        ATCIns = Allele.apply("ATC");

        NoCall = Allele.apply(".");
    }

    @Test
    public void testCreatingSNPAlleles() {
        Assert.assertTrue(A.isNonReference());
        Assert.assertFalse(A.isReference());
        Assert.assertTrue(A.basesMatch("A"));
        Assert.assertEquals(A.length(), 1);

        Assert.assertTrue(ARef.isReference());
        Assert.assertFalse(ARef.isNonReference());
        Assert.assertTrue(ARef.basesMatch("A"));
        Assert.assertFalse(ARef.basesMatch("T"));

        Assert.assertTrue(T.isNonReference());
        Assert.assertFalse(T.isReference());
        Assert.assertTrue(T.basesMatch("T"));
        Assert.assertFalse(T.basesMatch("A"));
    }

    @Test
    public void testCreatingNoCallAlleles() {
        Assert.assertTrue(NoCall.isNonReference());
        Assert.assertFalse(NoCall.isReference());
        Assert.assertFalse(NoCall.basesMatch("."));
        Assert.assertEquals(NoCall.length(), 0);
        Assert.assertTrue(NoCall.isNoCall());
        Assert.assertFalse(NoCall.isCalled());
    }


    @Test
    public void testCreatingIndelAlleles() {
        Assert.assertEquals(ATIns.length(), 2);
        Assert.assertEquals(ATCIns.length(), 3);
        Assert.assertEquals(ATIns.getBases(), "AT".getBytes());
        Assert.assertEquals(ATCIns.getBases(), "ATC".getBytes());
    }


    @Test
    public void testConstructors1() {
        Allele a1 = Allele.apply("A");
        Allele a2 = Allele.apply("A".getBytes());
        Allele a3 = Allele.apply("A");
        Allele a4 = Allele.apply("A", true);

        Assert.assertTrue(a1.equals(a2));
        Assert.assertTrue(a1.equals(a3));
        Assert.assertFalse(a1.equals(a4));
    }

    @Test
    public void testInsConstructors() {
        Allele a1 = Allele.apply("AC");
        Allele a2 = Allele.apply("AC".getBytes());
        Allele a3 = Allele.apply("AC");
        Allele a4 = Allele.apply("AC", true);

        Assert.assertTrue(a1.equals(a2));
        Assert.assertTrue(a1.equals(a3));
        Assert.assertFalse(a1.equals(a4));
    }

    @Test
    public void testEquals() {
        Assert.assertTrue(ARef.basesMatch(A));
        Assert.assertFalse(ARef.equals(A));
        Assert.assertFalse(ARef.equals(ATIns));
        Assert.assertFalse(ARef.equals(ATCIns));

        Assert.assertTrue(T.basesMatch(T));
        Assert.assertFalse(T.basesMatch(A));
        Assert.assertFalse(T.equals(A));

        Assert.assertTrue(ATIns.equals(ATIns));
        Assert.assertFalse(ATIns.equals(ATCIns));
        Assert.assertTrue(ATIns.basesMatch("AT"));
        Assert.assertFalse(ATIns.basesMatch("A"));
        Assert.assertFalse(ATIns.basesMatch("ATC"));

        Assert.assertTrue(ATIns.basesMatch("AT"));
        Assert.assertFalse(ATIns.basesMatch("ATC"));
    }

    @Test (expectedExceptions = IllegalArgumentException.class)
    public void testBadConstructorArgs1() {
        byte[] foo = null;
        Allele.apply(foo);
    }

    @Test (expectedExceptions = IllegalArgumentException.class)
    public void testBadConstructorArgs2() {
        Allele.apply("x");
    }

    @Test (expectedExceptions = IllegalArgumentException.class)
    public void testBadConstructorArgs3() {
        Allele.apply("--");
    }

    @Test (expectedExceptions = IllegalArgumentException.class)
    public void testBadConstructorArgs4() {
        Allele.apply("-A");
    }

    @Test (expectedExceptions = IllegalArgumentException.class)
    public void testBadConstructorArgs5() {
        Allele.apply("A A");
    }

    @Test
    public void testExtend() {
        Assert.assertEquals("AT", Allele.extend(Allele.apply("A"), "T".getBytes()).toString());
        Assert.assertEquals("ATA", Allele.extend(Allele.apply("A"), "TA".getBytes()).toString());
        Assert.assertEquals("A", Allele.extend(Allele.NO_CALL(), "A".getBytes()).toString());
        Assert.assertEquals("ATCGA", Allele.extend(Allele.apply("AT"), "CGA".getBytes()).toString());
        Assert.assertEquals("ATCGA", Allele.extend(Allele.apply("ATC"), "GA".getBytes()).toString());
    }
}

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



/**
 * Immutable representation of an allele
 *
 * Types of alleles:
 *
 * Ref: a t C g a // C is the reference base
 *
 *    : a t G g a // C base is a G in some individuals
 *
 *    : a t - g a // C base is deleted w.r.t. the reference
 *
 *    : a t CAg a // A base is inserted w.r.t. the reference sequence
 *
 * In these cases, where are the alleles?
 *
 * SNP polymorphism of C/G  -> { C , G } -> C is the reference allele
 * 1 base deletion of C     -> { tC , t } -> C is the reference allele and we include the preceding reference base (null alleles are not allowed)
 * 1 base insertion of A    -> { C ; CA } -> C is the reference allele (because null alleles are not allowed)
 *
 * Suppose I see a the following in the population:
 *
 * Ref: a t C g a // C is the reference base
 *    : a t G g a // C base is a G in some individuals
 *    : a t - g a // C base is deleted w.r.t. the reference
 *
 * How do I represent this?  There are three segregating alleles:
 *
 *  { C , G , - }
 *
 *  and these are represented as:
 *
 *  { tC, tG, t }
 *
 * Now suppose I have this more complex example:
 *
 * Ref: a t C g a // C is the reference base
 *    : a t - g a
 *    : a t - - a
 *    : a t CAg a
 *
 * There are actually four segregating alleles:
 *
 *   { Cg , -g, --, and CAg } over bases 2-4
 *
 *   represented as:
 *
 *   { tCg, tg, t, tCAg }
 *
 * Critically, it should be possible to apply an allele to a reference sequence to create the
 * correct haplotype sequence:
 *
 * Allele + reference => haplotype
 *
 * For convenience, we are going to create Alleles where the GenomeLoc of the allele is stored outside of the
 * Allele object itself.  So there's an idea of an A/C polymorphism independent of it's surrounding context.
 *
 * Given list of alleles it's possible to determine the "type" of the variation
 *
 *      A / C @ loc => SNP
 *      - / A => INDEL
 *
 * If you know where allele is the reference, you can determine whether the variant is an insertion or deletion.
 *
 * Alelle also supports is concept of a NO_CALL allele.  This Allele represents a haplotype that couldn't be
 * determined. This is usually represented by a '.' allele.
 *
 * Note that Alleles store all bases as bytes, in **UPPER CASE**.  So 'atc' == 'ATC' from the perspective of an
 * Allele.

 * @author ebanks, depristo
 */

/**
 * Ported to scala by Wim Spee on 2/20/14.
 */


package org.broadinstitute.variant.variantcontext

import java.util.Arrays

class Allele(private val _bases : Array[Byte],  val isReference : Boolean) extends Ordered[Allele] {


  //auxillary constructor, not sure if needed because of the constructors in the companion object
  protected def  this( bases : String,  isRef: Boolean)  {
    this(bases.getBytes(), isRef);
  }


  // null alleles are no longer allowed
  if ( Allele.wouldBeNullAllele(_bases) )      { throw new IllegalArgumentException("Null alleles are not supported")                                }
  if (! Allele.acceptableAlleleBases(_bases) ) { throw new IllegalArgumentException("Unexpected base in allele bases \'" + new String(_bases)+"\'");  }

  // ---------------------------------------------------------------------------------------------------------
  //
  // Allele members
  //
  // ---------------------------------------------------------------------------------------------------------

  // Returns true if this Allele is not the reference allele
  val isNonReference = !isReference

  //if this is the NO_CALL allele
  val isNoCall =  Allele.wouldBeNoCallAllele(_bases)
  val isCalled  = !isNoCall                           // Returns true if this is not the NO_CALL allele

  val bases = if(isNoCall){Allele.EMPTY_ALLELE_BASES}else{ _bases }

  // if this Allele is symbolic (i.e. no well-defined base sequence)
  val isSymbolic = Allele.wouldBeSymbolicAllele(bases)

  if(isNoCall && isReference)   { throw new IllegalArgumentException("Cannot tag a NoCall allele as the reference allele")  }
  if(isSymbolic && isReference) { throw new IllegalArgumentException("Cannot tag a symbolic allele as the reference allele")}

  //
  //    /**
  //     * Creates a new allele based on the provided one.  Ref state will be copied unless ignoreRefState is true
  //     * (in which case the returned allele will be non-Ref).
  //     *
  //     * This method is efficient because it can skip the validation of the bases (since the original allele was already validated)
  //     *
  //     * @param allele  the allele from which to copy the bases
  //     * @param ignoreRefState  should we ignore the reference state of the input allele and use the default ref state?
  //     */
  //    def Allele( allele : Allele,  ignoreRefState : Boolean) {
  //        this.bases = allele.bases;
  //        this.isRef = ignoreRefState ? false : allele.isRef;
  //        this.isNoCall = allele.isNoCall;
  //        this.isSymbolic = allele.isSymbolic;
  //    }



  // Returns a nice string representation of this object
  override def toString() : String =
  {
    val part1 = if( isNoCall ){ Allele.NO_CALL_STRING}else{getDisplayString()}
    val part2 = if(isReference ){"*"}else{""}
    part1 + part2
    //return ( isNoCall() ? NO_CALL_STRING : getDisplayString() ) + (isReference() ? "*" : "");
  }

  /**
   * Return the DNA bases segregating in this allele.  Note this isn't reference polarized,
   * so the Null allele is represented by a vector of length 0
   *
   * @return the segregating bases
   */
  def getBases() = { if(isSymbolic){Allele.EMPTY_ALLELE_BASES}else{bases} }
  // return isSymbolic ? EMPTY_ALLELE_BASES : bases;

  /**
   * Return the DNA bases segregating in this allele in String format.
   * This is useful, because toString() adds a '*' to reference alleles and getBases() returns garbage when you call toString() on it.
   *
   * @return the segregating bases
   */
  def getBaseString() =  { if(isNoCall ){ Allele.NO_CALL_STRING} else{ new String(getBases())} }
  //return isNoCall() ? NO_CALL_STRING : new String(getBases());
  /**
   * Return the printed representation of this allele.
   * Same as getBaseString(), except for symbolic alleles.
   * For symbolic alleles, the base string is empty while the display string contains <TAG>.
   *
   * @return the allele string representation
   */
  def getDisplayString() : String = { new String(bases); }

  /**
   * Same as #getDisplayString() but returns the result as byte[].
   *
   * Slightly faster then getDisplayString()
   *
   * @return the allele string representation
   */
  def getDisplayBases() = { bases; }

  /**
   * @param other  the other allele
   *
   * @return true if these alleles are equal
   */
  override def equals( other: Any) : Boolean= {
    if(other.isInstanceOf[Allele]){equals(other.asInstanceOf[Allele], false)}else{false}
  }

  /**
   * @return hash code
   */
  override def hashCode() = {
    var hash = 1;
    var limit = bases.length
    var counter = 0
    while(counter < limit)
    {
      hash += (counter+1) * bases(counter);
      counter +=1
    }
    hash
  }

  /**
   * Returns true if this and other are equal.  If ignoreRefState is true, then doesn't require both alleles has the
   * same ref tag
   *
   * @param other            allele to compare to
   * @param ignoreRefState   if true, ignore ref state in comparison
   * @return true if this and other are equal
   */
  def equals( other : Allele,  ignoreRefState : Boolean) = {
    (this eq other) || (isReference == other.isReference || ignoreRefState) && isNoCall == other.isNoCall && (bases == other.bases || Arrays.equals(bases, other.bases));
  }

  /**
   * @param test  bases to test against
   *
   * @return  true if this Allele contains the same bases as test, regardless of its reference status; handles Null and NO_CALL alleles
   */
  def basesMatch( test : Array[Byte]): Boolean  = {  !isSymbolic && (bases == test || Arrays.equals(bases, test)); }


  /**
   * @param test  bases to test against
   *
   * @return  true if this Allele contains the same bases as test, regardless of its reference status; handles Null and NO_CALL alleles
   */
  def basesMatch( test : String) : Boolean = {  basesMatch(test.toUpperCase().getBytes()); }

  /**
   * @param test  allele to test against
   *
   * @return  true if this Allele contains the same bases as test, regardless of its reference status; handles Null and NO_CALL alleles
   */
  def basesMatch( test: Allele) : Boolean =  {  basesMatch(test.getBases()); }

  /**
   * @return the length of this allele.  Null and NO_CALL alleles have 0 length.
   */
  def length() = {if(isSymbolic){0}else{bases.length} }

  override def compare( other : Allele) : Int = {
    if ( isReference  && other.isNonReference  )
      -1;
    else if ( isNonReference && other.isReference )
      1;
    else
      getBaseString().compareTo(other.getBaseString()); // todo -- potential performance issue
  }


}

object Allele
{
  private val EMPTY_ALLELE_BASES = new Array[Byte](0);

  /** A generic static NO_CALL allele for use */
  val NO_CALL_STRING : String = "."

  private val REF_A = new Allele("A", true)
  private val ALT_A = new Allele("A", false)
  private val REF_C = new Allele("C", true)
  private val ALT_C = new Allele("C", false)
  private val REF_G = new Allele("G", true)
  private val ALT_G = new Allele("G", false)
  private val REF_T = new Allele("T", true)
  private val ALT_T = new Allele("T", false)
  private val REF_N = new Allele("N", true)
  private val ALT_N = new Allele("N", false)
  val NO_CALL = new Allele(NO_CALL_STRING, false)


  // ---------------------------------------------------------------------------------------------------------
  //
  // creation routines
  //
  // ---------------------------------------------------------------------------------------------------------

  /**
   * Create a new Allele that includes bases and if tagged as the reference allele if isRef == true.  If bases
   * == '-', a Null allele is created.  If bases ==  '.', a no call Allele is created.
   *
   * @param bases the DNA sequence of this variation, '-', of '.'
   * @param isRef should we make this a reference allele?
   * @throws IllegalArgumentException if bases contains illegal characters or is otherwise malformated
   */
  def apply( bases  : Array[Byte],  isRef : Boolean) : Allele  = {
    if ( bases == null )
    {
      throw new IllegalArgumentException("create: the Allele base string cannot be null; use new Allele() or new Allele(\"\") to create a Null allele");
    }


    if ( bases.length == 1 ) {
      // optimization to return a static constant Allele for each single base object
      val returnAllele : Allele = bases(0) match
      {
        case '.' =>  if ( isRef ) throw new IllegalArgumentException("Cannot tag a NoCall allele as the reference allele");  NO_CALL;
        case 'A'| 'a' => (if(isRef) {REF_A}else{ALT_A})
        case 'C'| 'c' => (if(isRef) {REF_C}else{ALT_C})
        case 'G'| 'g' => (if(isRef) {REF_G}else{ALT_G})
        case 'T'| 't' => (if(isRef) {REF_T}else{ALT_T})
        case 'N'| 'n' => (if(isRef) {REF_N}else{ALT_N})
        case _ =>  throw new IllegalArgumentException("Illegal base [" + bases(0).asInstanceOf[Char] + "] seen in the allele");
      }
      returnAllele
    } else {
      new Allele(bases, isRef);
    }
  }

  def apply( base : Byte,  isRef : Boolean) : Allele =  { apply( Array[Byte](base), isRef)  }
  def apply( base : Byte) : Allele =                    { apply( base, false )              }

  /**
   * @see Allele(byte[], boolean)
   *
   * @param bases  bases representing an allele
   * @param isRef  is this the reference allele?
   */
  def apply( bases : String,  isRef : Boolean): Allele ={ apply(bases.getBytes(), isRef)    }


  /**
   * Creates a non-Ref allele.  @see Allele(byte[], boolean) for full information
   *
   * @param bases  bases representing an allele
   */
  def apply( bases : String) : Allele =                  { apply(bases, false)               }

  /**
   * Creates a non-Ref allele.  @see Allele(byte[], boolean) for full information
   *
   * @param bases  bases representing an allele
   */
  def apply( bases : Array[Byte]) : Allele =             { apply(bases, false)                }

//  /**
//   * Creates a new allele based on the provided one.  Ref state will be copied unless ignoreRefState is true
//   * (in which case the returned allele will be non-Ref).
//   *
//   * This method is efficient because it can skip the validation of the bases (since the original allele was already validated)
//   *
//   * @param allele  the allele from which to copy the bases
//   * @param ignoreRefState  should we ignore the reference state of the input allele and use the default ref state?
   //*/
  //    def apply(  allele : Allele,   ignoreRefState : Boolean) : Allele =  {
  //         new Allele(allele, ignoreRefState);
  //    }
  def extend( left : Allele,  right : Array[Byte]): Allele =  {
    if ( left.isSymbolic )
      throw new IllegalArgumentException("Cannot extend a symbolic allele");
    val bases = new Array[Byte](left.length() + right.length);
    System.arraycopy(left.getBases(), 0, bases, 0, left.length());
    System.arraycopy(right, 0, bases, left.length(), right.length);

    apply(bases, left.isReference );
  }

  /**
   * @param bases  bases representing an allele
   * @return true if the bases represent the null allele
   */
  def wouldBeNullAllele( bases : Array[Byte]) = {
    (bases.length == 1 && bases(0) == '-') || bases.length == 0;
  }

  /**
   * @param bases  bases representing an allele
   * @return true if the bases represent the NO_CALL allele
   */
  def wouldBeNoCallAllele( bases : Array[Byte]) : Boolean = {
    bases.length == 1 && bases(0) == '.';
  }

  /**
   * @param bases  bases representing an allele
   * @return true if the bases represent a symbolic allele
   */
  def wouldBeSymbolicAllele(bases : Array[Byte]) = {
    if ( bases.length <= 2 )
      false;
    else {
      val strBases = new String(bases);
      (bases(0) == '<' && bases(bases.length-1) == '>') || (strBases.contains("[") || strBases.contains("]"));
    }
  }

  /**
   * @param bases  bases representing an allele
   * @return true if the bases represent the well formatted allele
   */
  def acceptableAlleleBases( bases : String) : Boolean= {
    acceptableAlleleBases(bases.getBytes(), true);
  }

  def acceptableAlleleBases( bases : String,  allowNsAsAcceptable : Boolean) : Boolean = {
    acceptableAlleleBases(bases.getBytes(), allowNsAsAcceptable);
  }

  /**
   * @param bases  bases representing an allele
   * @return true if the bases represent the well formatted allele
   */
  def acceptableAlleleBases( bases : Array[Byte]) : Boolean =  {
    val acceptable = acceptableAlleleBases(bases, true); // default: N bases are acceptable
    acceptable
  }

  def acceptableAlleleBases( bases : Array[Byte],  allowNsAsAcceptable : Boolean) : Boolean = {
    if ( wouldBeNullAllele(bases) ){ false;}
    else
    {
      if ( wouldBeNoCallAllele(bases) || wouldBeSymbolicAllele(bases) ){true}
      else
      {
        var acceptable = true
        var counter = 0
        val limit = bases.length
        while(counter < limit)
        {
          bases(counter) match{
            case 'A'| 'C'| 'G'| 'T'| 'a'| 'c'| 'g'|'t' => {}
            case 'N' | 'n' => if(allowNsAsAcceptable){}else{acceptable = false}
            case _ => acceptable = false
          }
          counter += 1
        }
        acceptable
      }
    }
  }


  def getMatchingAllele(allAlleles : Array[Allele] , alleleBases : Array[Byte]) : Allele =  {
    if ( wouldBeNoCallAllele(alleleBases) ) { NO_CALL }
    else
    {
      var counter = 0
      val limit = allAlleles.size
      val iterator = allAlleles.iterator
      var allele : Allele = null
      while(counter < limit)
      {
        val a = iterator.next
        if ( a.basesMatch(alleleBases) ) {  allele =  a; }
        counter += 1
      }
      allele
    }
  }



  def oneIsPrefixOfOther( a1 : Allele, a2 : Allele)=  {
    if ( a2.length() >= a1.length() )
      firstIsPrefixOfSecond(a1, a2);
    else
      firstIsPrefixOfSecond(a2, a1);
  }

  def firstIsPrefixOfSecond( a1 : Allele, a2 : Allele)=  {
    val a1String = a1.getBaseString();
    a2.getBaseString().substring(0, a1String.length()).equals(a1String);
  }


}


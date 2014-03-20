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
 * Ported to Scala by Wim Spee on 2/20/14.
 */


package org.broadinstitute.variant.variantcontext

import com.google.java.contract._
import org.broadinstitute.variant.vcf.VCFConstants



import scala.collection._


/**
 * This class encompasses all the basic information about a genotype.  It is immutable.
 *
 * @author Mark DePristo
 */
//@Invariant({
//        "getAlleles() != null",
//        "getSampleName() != null",
//        "getPloidy() >= 0",
//        "! hasForbiddenKey(getExtendedAttributes())"})

abstract class Genotype(private val sampleName: String, private val filters: String) extends Comparable[Genotype] {

  private lazy val gType : GenotypeType = determineType()

  /**
   * @return the alleles for this genotype.  Cannot be null.  May be empty
   */
  @Ensures(Array[String]("result != null"))
  def getAlleles() : Array[Allele]

  /**
   * Returns how many times allele appears in this genotype object?
   *
   * @param allele
   * @return a value >= 0 indicating how many times the allele occurred in this sample's genotype
   */
  @Requires(Array[String]("allele != null"))
  @Ensures(Array[String]("result >= 0"))
  def countAllele(  allele : Allele) : Int = { getAlleles().count(_.equals(allele)) }

  /**
   * Get the ith allele in this genotype
   *
   * @param i the ith allele, must be < the ploidy, starting with 0
   * @return the allele at position i, which cannot be null
   */
  @Requires(Array[String]("i >=0 && i < getPloidy()","getType() != GenotypeType.UNAVAILABLE" ))
  @Ensures(Array[String]("result != null"))
  def getAllele( i : Int) : Allele

  /**
   * Are the alleles phased w.r.t. the global phasing system?
   *
   * @return true if yes
   */
  def isPhased() : Boolean

  /**
   * What is the ploidy of this sample?
   *
   * @return the ploidy of this genotype.  0 if the site is no-called.
   */
  @Ensures(Array[String]("result >= 0"))
  def getPloidy() : Int = { getAlleles().size;  }

  /**
   * @return the sequencing depth of this sample, or -1 if this value is missing
   */
  @Ensures(Array[String]("result >= -1"))
  def getDP() : Int

  /**
   * @return the count of reads, one for each allele in the surrounding Variant context,
   *      matching the corresponding allele, or null if this value is missing.  MUST
   *      NOT BE MODIFIED!
   */
  def getAD() : Array[Int]

  /**
   * Returns the name associated with this sample.
   *
   * @return a non-null String
   */
  @Ensures(Array[String]("result != null"))
  def getSampleName(): String =  {  sampleName;  }

  /**
   * Returns a phred-scaled quality score, or -1 if none is available
   * @return
   */
  @Ensures(Array[String]("result >= -1"))
  def getGQ() : Int

  /**
   * Does the PL field have a value?
   * @return true if there's a PL field value
   */
  @Ensures(Array[String]("(result == false && getPL() == null) || (result == true && getPL() != null)"))
  def hasPL() :Boolean = {  getPL() != null; }

  /**
   * Does the AD field have a value?
   * @return true if there's a AD field value
   */
  @Ensures(Array[String]("(result == false && getAD() == null) || (result == true && getAD() != null)"))
  def hasAD() : Boolean = {  getAD() != null;  }

  /**
   * Does the GQ field have a value?
   * @return true if there's a GQ field value
   */
  @Ensures(Array[String]("(result == false && getGQ() == -1) || (result == true && getGQ() >= 0)"))
  def hasGQ() : Boolean = { getGQ() != -1; }

  /**
   * Does the DP field have a value?
   * @return true if there's a DP field value
   */
  @Ensures(Array[String]("(result == false && getDP() == -1) || (result == true && getDP() >= 0)"))
  def hasDP() : Boolean =  { getDP() != -1; }

  // ---------------------------------------------------------------------------------------------------------
  //
  // The type of this genotype
  //
  // ---------------------------------------------------------------------------------------------------------

  /**
   * @return the high-level type of this sample's genotype
   */
  @Ensures(Array[String]("type != null", "result != null"))
  def getType() : GenotypeType = { gType }

  /**
   * Internal code to determine the type of the genotype from the alleles vector
   * @return the type
   */
  @Requires(Array[String]("type == null")) // we should never call if already calculated
  protected def determineType() : GenotypeType = {

    val alleles = getAlleles();
    if ( alleles.isEmpty ){  GenotypeType.UNAVAILABLE;}
    else
    {
      var sawNoCall = false
      var sawMultipleAlleles = false;
      var observedAllele : Allele= null;

      var counter = 0
      val limit = alleles.size

      while(counter < limit)
      {
        val allele = alleles(counter)
        if ( allele.isNoCall ){   sawNoCall = true;}
        else if ( observedAllele == null ){  observedAllele = allele;}
        else if ( !allele.equals(observedAllele) ){ sawMultipleAlleles = true;}
        counter +=1
      }

      if ( sawNoCall ) {
        if ( observedAllele == null ){GenotypeType.NO_CALL} else{GenotypeType.MIXED}

      }
      else
      {
        if ( observedAllele == null ){ throw new IllegalStateException("BUG: there are no alleles present in this genotype but the alleles list is not null");}
        else
        {
          if( sawMultipleAlleles){ GenotypeType.HET}
          else
          {
            if(observedAllele.isReference){ GenotypeType.HOM_REF}  else{ GenotypeType.HOM_VAR}
          }
        }

      }



    }

    //        // TODO -- this code is slow and could be optimized for the diploid case
    //        final List<Allele> alleles = getAlleles();
    //        if ( alleles.isEmpty() )
    //            return GenotypeType.UNAVAILABLE;
    //
    //        boolean sawNoCall = false, sawMultipleAlleles = false;
    //        Allele observedAllele = null;
    //
    //        for ( final Allele allele : alleles ) {
    //            if ( allele.isNoCall() )
    //                sawNoCall = true;
    //            else if ( observedAllele == null )
    //                observedAllele = allele;
    //            else if ( !allele.equals(observedAllele) )
    //                sawMultipleAlleles = true;
    //        }
    //
    //        if ( sawNoCall ) {
    //            if ( observedAllele == null )
    //                return GenotypeType.NO_CALL;
    //            return GenotypeType.MIXED;
    //        }
    //
    //        if ( observedAllele == null )
    //            throw new IllegalStateException("BUG: there are no alleles present in this genotype but the alleles list is not null");
    //
    //        return sawMultipleAlleles ? GenotypeType.HET : observedAllele.isReference() ? GenotypeType.HOM_REF : GenotypeType.HOM_VAR;
  }

  /**
   * @return true if all observed alleles are the same (regardless of whether they are ref or alt); if any alleles are no-calls, this method will return false.
   */
  def isHom() : Boolean =  { isHomRef() || isHomVar(); }

  /**
   * @return true if all observed alleles are ref; if any alleles are no-calls, this method will return false.
   */
  def isHomRef() : Boolean = { gType == GenotypeType.HOM_REF; }

  /**
   * @return true if all observed alleles are alt; if any alleles are no-calls, this method will return false.
   */
  def isHomVar(): Boolean = {  gType == GenotypeType.HOM_VAR; }

  /**
   * @return true if we're het (observed alleles differ); if the ploidy is less than 2 or if any alleles are no-calls, this method will return false.
   */
  def isHet() : Boolean =  { gType == GenotypeType.HET; }

  /**
   * @return true if this genotype is not actually a genotype but a "no call" (e.g. './.' in VCF); if any alleles are not no-calls (even if some are), this method will return false.
   */
  def isNoCall() : Boolean = { gType == GenotypeType.NO_CALL; }

  /**
   * @return true if this genotype is comprised of any alleles that are not no-calls (even if some are).
   */
  def isCalled() : Boolean =  { gType != GenotypeType.NO_CALL && gType != GenotypeType.UNAVAILABLE; }

  /**
   * @return true if this genotype is comprised of both calls and no-calls.
   */
  def isMixed() : Boolean = { gType == GenotypeType.MIXED; }

  /**
   * @return true if the type of this genotype is set.
   */
  def isAvailable() : Boolean =  { gType != GenotypeType.UNAVAILABLE; }

  // ------------------------------------------------------------------------------
  //
  // methods for getting genotype likelihoods for a genotype object, if present
  //
  // ------------------------------------------------------------------------------

  /**
   * @return Returns true if this Genotype has PL field values
   */
  @Ensures(Array[String]("(result && getLikelihoods() != null) || (! result && getLikelihoods() == null)"))
  def  hasLikelihoods() : Boolean =  { getPL() != null;  }

  /**
   * Convenience function that returns a string representation of the PL field of this
   * genotype, or . if none is available.
   *
   * @return a non-null String representation for the PL of this sample
   */
  @Ensures(Array[String]("result != null"))
  def getLikelihoodsString() : String = {
    if(hasLikelihoods()){ getLikelihoods().toString()} else{ VCFConstants.MISSING_VALUE_v4};
  }

  /**
   * Returns the GenotypesLikelihoods data associated with this Genotype, or null if missing
   * @return null or a GenotypesLikelihood object for this sample's PL field
   */
  @Ensures(Array[String]("(hasLikelihoods() && result != null) || (! hasLikelihoods() && result == null)"))
  def getLikelihoods() : GenotypeLikelihoods =  {
    if(hasLikelihoods()){ GenotypeLikelihoods.fromPLs(getPL())} else{ null;}
  }

  /**
   * Are all likelihoods for this sample non-informative?
   *
   * Returns true if all PLs are 0 => 0,0,0 => true
   * 0,0,0,0,0,0 => true
   * 0,10,100 => false
   *
   * @return true if all samples PLs are equal and == 0
   */
  def isNonInformative() : Boolean = {
    if ( getPL() == null ){ return true;}
    else {
      val plArray = getPL()
      var counter = 0
      val limit = plArray.size
      var isNonInformativeBool = true

      while(counter < limit)
      {
        val PL = plArray(counter)
        if ( PL != 0 ){isNonInformativeBool = false}
        counter +=1
      }
      isNonInformativeBool
    }
  }

  /**
   * Unsafe low-level accessor the PL field itself, may be null.
   *
   * @return a pointer to the underlying PL data.  MUST NOT BE MODIFIED!
   */
  def getPL() : Array[Int]

  // ---------------------------------------------------------------------------------------------------------
  //
  // Many different string representations
  //
  // ---------------------------------------------------------------------------------------------------------

  /**
   * Return a VCF-like string representation for the alleles of this genotype.
   *
   * Does not append the reference * marker on the alleles.
   *
   * @return a string representing the genotypes, or null if the type is unavailable.
   */
  @Ensures(Array[String]("result != null || ! isAvailable()"))
  def getGenotypeString() : String =  {
    getGenotypeString(true);
  }

  /**
   * Return a VCF-like string representation for the alleles of this genotype.
   *
   * If ignoreRefState is true, will not append the reference * marker on the alleles.
   *
   * @return a string representing the genotypes, or null if the type is unavailable.
   */


  @Ensures(Array[String]("result != null || ! isAvailable()"))
  def getGenotypeString( ignoreRefState : Boolean) : String = {
    if ( getPloidy() == 0 ){ "NA";}
    // Notes:
    // 1. Make sure to use the appropriate separator depending on whether the genotype is phased
    // 2. If ignoreRefState is true, then we want just the bases of the Alleles (ignoring the '*' indicating a ref Allele)
    // 3. So that everything is deterministic with regards to integration tests, we sort Alleles (when the genotype isn't phased, of course)
    else
    {
      val SEPARATOR = if(isPhased()) { Genotype.PHASED_ALLELE_SEPARATOR} else{Genotype.UNPHASED_ALLELE_SEPARATOR}
      if(ignoreRefState){getAlleleStrings().mkString(SEPARATOR)}
      else
      {
        if(isPhased()){getAlleles().mkString(SEPARATOR)}
        else{Genotype.sortAllelesToNewArray(getAlleles()).mkString(SEPARATOR)}
      }
    }
  }

  /**
   * Utility that returns a list of allele strings corresponding to the alleles in this sample
   * @return
   */
  protected def getAlleleStrings() : Array[String] = { getAlleles().map(_.getBaseString()) }

  override def toString() : String = {
    String.format("[%s %s%s%s%s%s%s%s]",
      getSampleName(),
      getGenotypeString(false),
      Genotype.toStringIfExists(VCFConstants.GENOTYPE_QUALITY_KEY, getGQ()),
      Genotype.toStringIfExists(VCFConstants.DEPTH_KEY, getDP()),
      Genotype.toStringIfExists(VCFConstants.GENOTYPE_ALLELE_DEPTHS, getAD()),
      Genotype.toStringIfExists(VCFConstants.GENOTYPE_PL_KEY, getPL()),
      Genotype.toStringIfExists(VCFConstants.GENOTYPE_FILTER_KEY, getFilters()),
      Genotype.sortedString(getExtendedAttributes()));
  }

  def toBriefString() : String = {
    String.format("%s:Q%d", getGenotypeString(false), new Integer(getGQ()));
  }

  // ---------------------------------------------------------------------------------------------------------
  //
  // Comparison operations
  //
  // ---------------------------------------------------------------------------------------------------------

  /**
   * comparable genotypes -> compareTo on the sample names
   * @param genotype
   * @return
   */
  override
  def compareTo(  genotype : Genotype) : Int =  {
    getSampleName().compareTo(genotype.getSampleName());
  }

  def  sameGenotype(  other : Genotype)  : Boolean = {
    sameGenotype(other, true);
  }

  def sameGenotype( other : Genotype,  ignorePhase : Boolean) = {
    if (getPloidy() != other.getPloidy()){ false;} // gotta have the same number of allele to be equal
    else
    {
      // By default, compare the elements in the lists of alleles, element-by-element
      val thisAlleles : Array[Allele] = this.getAlleles();
      val otherAlleles : Array[Allele] = other.getAlleles();

      if (ignorePhase) { // do not care about order, only identity of Alleles
      val thisAllelesSorted = Genotype.sortAllelesToNewArray(thisAlleles)
        val otherAllelesSorted = Genotype.sortAllelesToNewArray(otherAlleles)
        thisAllelesSorted.sameElements(otherAllelesSorted)
      }
      else{  thisAlleles.sameElements(otherAlleles);}
    }
  }

  // ---------------------------------------------------------------------------------------------------------
  //
  // get routines for extended attributes
  //
  // ---------------------------------------------------------------------------------------------------------

  /**
   * Returns the extended attributes for this object
   * @return is never null, but is often isEmpty()
   */
  @Ensures(Array[String]("result != null", "! hasForbiddenKey(result)"))
  def getExtendedAttributes() : immutable.Map[String, Any]

  /**
   * Is key associated with a value (even a null one) in the extended attributes?
   *
   * Note this will not return true for the inline attributes DP, GQ, AD, or PL
   *
   * @param key a non-null string key to check for an association
   * @return true if key has a value in the extendedAttributes
   */
  @Requires(Array[String]("key != null", "! isForbiddenKey(key)"))
  def hasExtendedAttribute(  key : String)  : Boolean = {
    getExtendedAttributes().contains(key);
  }

  /**
   * Get the extended attribute value associated with key, if possible
   *
   * @param key a non-null string key to fetch a value for
   * @param defaultValue the value to return if key isn't in the extended attributes
   * @return a value (potentially) null associated with key, or defaultValue if no association exists
   */
  @Requires(Array[String]("key != null","! isForbiddenKey(key)" ))
  @Ensures(Array[String]("hasExtendedAttribute(key) || result == defaultValue"))
  def getExtendedAttribute( key : String,  defaultValue: Any) : Any = {
    if(hasExtendedAttribute(key)) { getExtendedAttributes()(key) } else{ defaultValue;}
  }

  /**
   * Same as #getExtendedAttribute with a null default
   *
   * @param key
   * @return
   */
  def getExtendedAttribute(key : String ) : Any = {
    getExtendedAttribute(key, null);
  }

  /**
   * Returns the filter string associated with this Genotype.
   *
   * @return If this result == null, then the genotype is considered PASSing filters
   *   If the result != null, then the genotype has failed filtering for the reason(s)
   *   specified in result.  To be reference compliant multiple filter field
   *   string values can be encoded with a ; separator.
   */
  def getFilters() = { filters  }

  /**
   * Is this genotype filtered or not?
   *
   * @return returns false if getFilters() == null
   */
  @Ensures(Array[String]("result != (getFilters() == null)"))
  def isFiltered() = { getFilters() != null;  }

  @Deprecated def hasLog10PError() =  { hasGQ(); }
  @Deprecated def getLog10PError() = {  getGQ() / -10.0; }
  @Deprecated def getPhredScaledQual() = { getGQ(); }

  @Deprecated
  def getAttributeAsString( key : String,  defaultValue : String) = {
    val x : Any = getExtendedAttribute(key);
    if ( x == null ){ defaultValue;}
    else
    {
      if ( x.isInstanceOf[String] ){ x.toString;}
      else(String.valueOf(x))
    }
  }

  @Deprecated
  def getAttributeAsInt( key : String, defaultValue : Int)  : Int = {
    val x : Any = getExtendedAttribute(key);
    if ( x == null || x == VCFConstants.MISSING_VALUE_v4 ) { defaultValue;}
    else
    {
      if ( x.isInstanceOf[Integer] ) { x.asInstanceOf[Integer]}
      else{ x.asInstanceOf[String].toInt;} // throws an exception if this isn't a string

    }

  }

  @Deprecated
  def getAttributeAsDouble( key : String,  defaultValue : Double) : Double = {
    val x : Any = getExtendedAttribute(key);
    if ( x == null ) { defaultValue;}
    else
    {
      if ( x.isInstanceOf[Double] ){ x.asInstanceOf[Double]}
      else{ x.asInstanceOf[String].toDouble} // throws an exception if this isn't a string
    }

  }

  /**
   * A totally generic getter, that allows you to specific keys that correspond
   * to even inline values (GQ, for example).  Can be very expensive.  Additionally,
   * all int[] are converted inline into List<Integer> for convenience.
   *
   * @param key
   * @return
   */
  def getAnyAttribute(key :  String ) : Any = {
    key match
    {
      case VCFConstants.GENOTYPE_KEY            => getAlleles()
      case VCFConstants.GENOTYPE_QUALITY_KEY    => getGQ()
      case VCFConstants.GENOTYPE_ALLELE_DEPTHS  => getAD()
      case VCFConstants.GENOTYPE_PL_KEY         => getPL()
      case VCFConstants.DEPTH_KEY               => getDP()
      case _                                    => getExtendedAttribute(key)
    }
  }

  def hasAnyAttribute(  key : String) : Boolean =  {
    if (key.equals(VCFConstants.GENOTYPE_KEY)) {
      return isAvailable();
    } else if (key.equals(VCFConstants.GENOTYPE_QUALITY_KEY)) {
      return hasGQ();
    } else if (key.equals(VCFConstants.GENOTYPE_ALLELE_DEPTHS)) {
      return hasAD();
    } else if (key.equals(VCFConstants.GENOTYPE_PL_KEY)) {
      return hasPL();
    } else if (key.equals(VCFConstants.DEPTH_KEY)) {
      return hasDP();
    } else {
      return hasExtendedAttribute(key);
    }
  }

  // TODO -- add getAttributesAsX interface here




}

object Genotype
{

  /**
   * A list of genotype field keys corresponding to values we
   * manage inline in the Genotype object.  They must not appear in the
   * extended attributes map
   */

  val PRIMARY_KEYS = Array(
    VCFConstants.GENOTYPE_FILTER_KEY,
    VCFConstants.GENOTYPE_KEY,
    VCFConstants.GENOTYPE_QUALITY_KEY,
    VCFConstants.DEPTH_KEY,
    VCFConstants.GENOTYPE_ALLELE_DEPTHS,
    VCFConstants.GENOTYPE_PL_KEY);

  val PHASED_ALLELE_SEPARATOR = "|";
  val UNPHASED_ALLELE_SEPARATOR = "/";

  /**
   * Does the attribute map have a mapping involving a forbidden key (i.e.,
   * one that's managed inline by this Genotypes object?
   *
   * @param attributes the extended attributes key
   * @return
   */
  protected def hasForbiddenKey( attributes : Map[String, Any] ) = {

    val forbiddenKeyCount = PRIMARY_KEYS.map(attributes.contains(_)).count(_ == true)
    val hasForbiddenKey = forbiddenKeyCount > 0
    hasForbiddenKey
  }

  protected def isForbiddenKey( key : String) : Boolean = { PRIMARY_KEYS.contains(key)}

  // ------------------------------------------------------------------------------
  //
  // private utilities
  //
  // ------------------------------------------------------------------------------



  /**
   * a utility method for generating sorted strings from a map key set.
   * @param map the map
   * @return a sting, enclosed in {}, with comma seperated key value pairs in order of the keys
   */
  @Requires(Array[String]("c != null"))
  protected def sortedString (map : Map[String, Any]) : String =  {

    // NOTE -- THIS IS COPIED FROM GATK UTILS TO ALLOW US TO KEEP A SEPARATION BETWEEN THE GATK AND VCF CODECS
    val keys : Array[String] = map.keySet.toArray.sorted
    val pairs = keys.map( key => key + "=" + map(key))

    if(pairs.isEmpty) { ""} else{ " {" + pairs.mkString(" ,") + "}"};


    //    val t = new ArrayList[T](c.keySet());
    //    Collections.sort(t);
    //
    //    val pairs = new ArrayList[String]();
    //    var counter = 0
    //    val limit = t.size()
    //    while(counter < limit)
    //    {
    //      val k = t.get(counter)
    //      pairs.add(k + "=" + c.get(k))
    //      counter += 1
    //    }
    //
    //    if(pairs.isEmpty()) { ""} else{ " {" + ParsingUtils.join(", ", pairs.toArray(new Array[String](pairs.size))) + "}"};

  }

  /**
   * Returns a display name for field name with value v if this isn't -1.  Otherwise returns ""
   * @param name of the field ("AD")
   * @param v the value of the field, or -1 if missing
   * @return a non-null string for display if the field is not missing
   */
  @Requires(Array[String]("name != null"))
  @Ensures(Array[String]("result != null"))
  protected def toStringIfExists( name : String,  v : Int) =  {  if( v == -1){ ""}else{ " " + name + " " + v}; }

  /**
   * Returns a display name for field name with String value v if this isn't null.  Otherwise returns ""
   * @param name of the field ("FT")
   * @param v the value of the field, or null if missing
   * @return a non-null string for display if the field is not missing
   */
  protected def toStringIfExists( name : String,  v : String) = { if(v == null){ ""} else{ " " + name + " " + v};}


  /**
   * Returns a display name for field name with values vs if this isn't null.  Otherwise returns ""
   * @param name of the field ("AD")
   * @param vs the value of the field, or null if missing
   * @return a non-null string for display if the field is not missing
   */
  @Requires(Array[String]("name != null"))
  @Ensures(Array[String]("result != null"))
  protected def toStringIfExists(name: String, vs: Array[Int]) = {
    if (vs == null) {
      ""
    }
    else {
      val b = new StringBuilder();
      b.append(" ").append(name).append(" ");
      var counter = 0
      val limit = vs.size
      while (counter < limit) {
        if (counter != 0) b.append(",");
        b.append(vs(counter));
        counter += 1

      }
      b.toString();
    }
  }

  def sortAllelesToNewArray(alleles : Array[Allele]) =
  {
    val sortedArray : Array[Allele] = new Array[Allele](alleles.size)
    alleles.copyToArray(sortedArray)
    sortedArray.sortWith(_ < _)
  }








}






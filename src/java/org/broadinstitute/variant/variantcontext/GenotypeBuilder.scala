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

package org.broadinstitute.variant.variantcontext

import com.google.java.contract.Ensures;
import com.google.java.contract.Invariant;
import com.google.java.contract.Requires;
import org.broadinstitute.variant.vcf.VCFConstants

import scala.collection._

/**
 * A builder class for genotypes
 *
 * Provides convenience setter methods for all of the Genotype field
 * values.  Setter methods can be used in any order, allowing you to
 * pass through states that wouldn't be allowed in the highly regulated
 * immutable Genotype class.
 *
 * All fields default to meaningful MISSING values.
 *
 * Call make() to actually create the corresponding Genotype object from
 * this builder.  Can be called multiple times to create independent copies,
 * or with intervening sets to conveniently make similar Genotypes with
 * slight modifications.
 *
 * @author Mark DePristo
 * @since 06/12
 */

/**
 * Created by Wim Spee on 2/23/14.
 */
@Invariant(Array[String]("alleles != null"))
/**
 * Create a empty builder.  Both a sampleName and alleles must be provided
 * before trying to make a Genotype from this builder.
 */
class GenotypeBuilder {

  private var sampleName : String = null;
  private var alleles : Array[Allele]= Array[Allele]()

  private var isPhased = false;
  private var GQ = -1;
  private var DP = -1;
  private var AD  : Array[Int]= null;
  private var PL : Array[Int] = null;
  private var extendedAttributes : immutable.HashMap[String, Any] = GenotypeBuilder.NO_ATTRIBUTES;
  private var filters : String = null;
  private var initialAttributeMapSize = 5;



  /**
   * Create a builder using sampleName.  Alleles must be provided
   * before trying to make a Genotype from this builder.
   * @param sampleName
   */
  def this( sampleName : String) {
    this()
    name(sampleName);
  }

  /**
   * Make a builder using sampleName and alleles for starting values
   * @param sampleName
   * @param allelesArgument
   */
  def this( sampleName : String,  allelesArgument :  Array[Allele]) {
    this()
    name(sampleName);
    alleles(allelesArgument);
  }

  /**
   * Create a new builder starting with the values in Genotype g
   * @param g
   */
  def this( g : Genotype) {
    this()
    copy(g);
  }

  /**
   * Copy all of the values for this builder from Genotype g
   * @param g
   * @return
   */
  def copy(  g : Genotype) : GenotypeBuilder =  {
    name(g.getSampleName());
    alleles(g.getAlleles());
    phased(g.isPhased());
    GQ(g.getGQ());
    DP(g.getDP());
    AD(g.getAD());
    PL(g.getPL());
    filter(g.getFilters());
    attributes(g.getExtendedAttributes());
    return this;
  }

  /**
   * Reset all of the builder attributes to their defaults.  After this
   * function you must provide sampleName and alleles before trying to
   * make more Genotypes.
   */
  def reset( keepSampleName : Boolean) {
    if ( ! keepSampleName ) sampleName = null;
    alleles = Array[Allele]();
    isPhased = false;
    GQ = -1;
    DP = -1;
    AD = null;
    PL = null;
    filters = null;
    extendedAttributes = null;
  }

  /**
   * Create a new Genotype object using the values set in this builder.
   *
   * After creation the values in this builder can be modified and more Genotypes
   * created, althrough the contents of array values like PL should never be modified
   * inline as they are not copied for efficiency reasons.
   *
   * @return a newly minted Genotype object with values provided from this builder
   */
  @Ensures(Array[String]("result != null"))
  def make() : Genotype =  {
    val ea = if(extendedAttributes == null){ GenotypeBuilder.NO_ATTRIBUTES} else{ extendedAttributes}
    new FastGenotype(sampleName, alleles, isPhased, GQ, DP, AD, PL, filters, ea);
  }

  /**
   * Set this genotype's name
   * @param sampleName
   * @return
   */
  @Requires(Array[String]("sampleName != null"))
  @Ensures(Array[String]("this.sampleName != null"))
  def name( sampleName : String) : GenotypeBuilder= {
    this.sampleName = sampleName;
    this;
  }

  /**
   * Set this genotype's alleles
   * @param alleles
   * @return
   */
  @Ensures(Array[String]("this.alleles != null"))
  def alleles( alleles :  Array[Allele]) : GenotypeBuilder =  {
    if ( alleles == null ){ this.alleles = Array[Allele]();}
    else{  this.alleles = alleles}
    this;
  }

  /**
   * Is this genotype phased?
   * @param phased
   * @return
   */
  def phased(  phased : Boolean) : GenotypeBuilder = {
    isPhased = phased;
    this;
  }

  @Requires(Array[String]("GQ >= -1"))
  @Ensures(Array[String]("this.GQ == GQ", "this.GQ >= -1"))
  def GQ( GQ : Int) : GenotypeBuilder = {
    this.GQ = GQ;
    return this;
  }

  /**
   * Adaptor interface from the pLog10Error system.
   *
   * Will be retired when
   *
   * @param pLog10Error
   * @return
   */
  @Deprecated
  def log10PError(  pLog10Error :Double) :GenotypeBuilder = {
    if ( pLog10Error == CommonInfo.NO_LOG10_PERROR )
    { GQ(-1);}
    else
    {GQ(Math.round(pLog10Error * -10).toInt);}
  }

  /**
   * This genotype has no GQ value
   * @return
   */
  def noGQ() :GenotypeBuilder =  { GQ = -1; this; }

  /**
   * This genotype has no AD value
   * @return
   */
  def noAD() : GenotypeBuilder = { AD = null; this; }

  /**
   * This genotype has no DP value
   * @return
   */
  def noDP() :GenotypeBuilder = { DP = -1;  this; }

  /**
   * This genotype has no PL value
   * @return
   */
  def noPL() : GenotypeBuilder =  { PL = null; this; }

  /**
   * This genotype has this DP value
   * @return
   */
  @Requires(Array[String]("DP >= -1"))
  @Ensures(Array[String]("this.DP == DP"))
  def DP(  DP :Int)  :  GenotypeBuilder = {
    this.DP = DP;
    this;
  }

  /**
   * This genotype has this AD value
   * @return
   */
  @Requires(Array[String]("AD == null || AD.length > 0"))
  @Ensures(Array[String]("this.AD == AD"))
  def AD( AD : Array[Int]) : GenotypeBuilder =  {
    this.AD = AD;
    this;
  }

  /**
   * This genotype has this PL value, as int[].  FAST
   * @return
   */
  @Requires(Array[String]("PL == null || PL.length > 0"))
  @Ensures(Array[String]("this.PL == PL"))
  def  PL( PL : Array[Int]) :GenotypeBuilder = {
    this.PL = PL;
    this;
  }

  /**
   * This genotype has this PL value, converted from double[]. SLOW
   * @return
   */
  @Requires(Array[String]("PL == null || PL.length > 0"))
  @Ensures(Array[String]({"this.PL == PL"}))
  def PL( GLs : Array[Double]) : GenotypeBuilder = {
    this.PL = GenotypeLikelihoods.fromLog10Likelihoods(GLs).getAsPLs();
    this;
  }

  /**
   * This genotype has these attributes.
   *
   * Cannot contain inline attributes (DP, AD, GQ, PL)
   * @return
   */
  @Requires(Array[String]("attributes != null"))
  @Ensures(Array[String]("attributes.isEmpty() || extendedAttributes != null"))
  def attributes(  attributes : immutable.Map[String, Any]) : GenotypeBuilder = {
    for((k,v) <- attributes){attribute(k,v)}

    this;
  }

  /**
   * This genotype has these attributes.
   *
   * Cannot contain inline attributes (DP, AD, GQ, PL)
   * @return
   */
  @Requires(Array[String]("attributes != null"))
  @Ensures(Array[String]("attributes.isEmpty() || extendedAttributes != null"))
  def attributesMapFromJava(  attributesMap : Map[String, Any]) : GenotypeBuilder = {
    attributes(attributesMap.toMap)
  }



  /**
   * Tells this builder to remove all extended attributes
   *
   * @return
   */
  def  noAttributes() : GenotypeBuilder = {
    this.extendedAttributes = null;
    this;
  }

  /**
   * This genotype has this attribute key / value pair.
   *
   * Cannot contain inline attributes (DP, AD, GQ, PL)
   * @return
   */
  @Requires(Array[String]({"key != null"}))
  @Ensures(Array[String]("extendedAttributes != null", "extendedAttributes.containsKey(key)"))
  def attribute(  key : String,   value : Any) : GenotypeBuilder=  {

    if(key == "GV")
    {
      val blaat = "blaat"
    }

    if(value.isInstanceOf[String] && value.asInstanceOf[String] == "Some(S1)")
    {
      val blaat = "blaat";
    }


    if ( extendedAttributes == null ){
      extendedAttributes = immutable.HashMap[String, Any]()}
    extendedAttributes = extendedAttributes + (key -> value)
    this;
  }

  /**
   * Tells this builder to make a Genotype object that has had filters applied,
   * which may be empty (passes) or have some value indicating the reasons
   * why it's been filtered.
   *
   * @param filters non-null list of filters.  empty list => PASS
   * @return this builder
   */
  @Requires(Array[String]("filters != null"))
  def filters(  filters : Array[String]) : GenotypeBuilder = {
    if ( filters.isEmpty )
    { filter(null);}
    else if ( filters.size == 1 )
    { filter(filters(0));}
    else
    {
      filter(filters.sorted.mkString(";"));
    }
  }

  /**
   * varargs version of #filters
   * @param filtersArgument
   * @return
   */
  @Requires(Array[String]("filters != null"))
  def filters( filtersArgument : String*) : GenotypeBuilder = {
    val size = filtersArgument.size
    val list = new Array[String](size)
    var counter = 0

    while(counter  < size)
    {
      list(counter) = (filtersArgument(counter))
      counter +=1
    }

    filters(list);
  }

  /**
   * Most efficient version of setting filters -- just set the filters string to filters
   *
   * @param filter if filters == null or filters.equals("PASS") => genotype is PASS
   * @return
   */
  def  filter( filter : String) : GenotypeBuilder =  {
    this.filters = if(VCFConstants.PASSES_FILTERS_v4.equals(filter)){ null}else{ filter}
    this;
  }

  /**
   * This genotype is unfiltered
   *
   * @return
   */
  def  unfiltered() : GenotypeBuilder =  { filter(null);
  }

  /**
   * Tell's this builder that we have at most these number of attributes
   * @return
   */
  def  maxAttributes( i : Int) : GenotypeBuilder = {
    initialAttributeMapSize = i;
    this;
  }
}

object GenotypeBuilder
{
  val HAPLOID_NO_CALL = Array[Allele](Allele.NO_CALL);
  val DIPLOID_NO_CALL = Array[Allele](Allele.NO_CALL, Allele.NO_CALL);

  val NO_ATTRIBUTES = new immutable.HashMap[String, Any]()


  // -----------------------------------------------------------------
  //
  // Factory methods
  //
  // -----------------------------------------------------------------

  def apply( sampleName : String,   alleles : Array[Allele])  : Genotype = {  new GenotypeBuilder(sampleName, alleles).make()  }



  def apply( sampleName : String, alleles  : Array[Allele], attributes : immutable.Map[String, Any]) : Genotype = {
    new GenotypeBuilder(sampleName, alleles).attributes(attributes).make();
  }

  def apply( sampleName : String,   alleles :  Array[Allele],  gls : Array[Double]) :Genotype  = {
    new GenotypeBuilder(sampleName, alleles).PL(gls).make();
  }

  /**
   * Create a new Genotype object for a sample that's missing from the VC (i.e., in
   * the output header).  Defaults to a diploid no call genotype ./.
   *
   * @param sampleName the name of this sample
   * @return an initialized Genotype with sampleName that's a diploid ./. no call genotype
   */
  def createMissing( sampleName : String, ploidy :Int) : Genotype =  {
    val builder = new GenotypeBuilder(sampleName);
    ploidy match {
      case 1 =>   builder.alleles(HAPLOID_NO_CALL);
      case 2 =>   builder.alleles(DIPLOID_NO_CALL);
      case _ => builder.alleles(Array.fill(ploidy){Allele.NO_CALL});
    }
    builder.make();
  }



}
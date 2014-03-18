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
 * Created by Wim Spee on 3/2/14.
 */


package org.broadinstitute.variant.variantcontext

import com.google.java.contract._
import org.broadinstitute.variant.vcf.VCFConstants
import scala.collection._
import scala.collection.JavaConversions._



/**
 * Builder class for VariantContext
 *
 * Some basic assumptions here:
 *
 * 1 -- data isn't protectively copied.  If you provide an attribute map to
 * the build, and modify it later, the builder will see this and so will any
 * resulting variant contexts.  It's best not to modify collections provided
 * to a builder.
 *
 * 2 -- the system uses the standard builder model, allowing the simple construction idiom:
 *
 *   builder.source("a").genotypes(gc).id("x").make() => VariantContext
 *
 * 3 -- The best way to copy a VariantContext is:
 *
 *   new VariantContextBuilder(vc).make() => a copy of VC
 *
 * 4 -- validation of arguments is done at the during the final make() call, so a
 * VariantContextBuilder can exist in an inconsistent state as long as those issues
 * are resolved before the call to make() is issued.
 *
 * @author depristo
 */

/**
 * Create an empty VariantContextBuilder where all values adopt their default values.  Note that
 * source, chr, start, stop, and alleles must eventually be filled in, or the resulting VariantContext
 * will throw an error.
 */
class VariantContextBuilder {
  // required fields
  private var fullyDecoded = false
  private var source : String = null
  private var contig : String = null
  private var start : Long = -1
  private var stop : Long = -1l
  private var alleleContext : AlleleContext = null

  // optional -> these are set to the appropriate default value
  private var ID = VCFConstants.EMPTY_ID_FIELD
  private var genotypeContext = GenotypesContext.NO_GENOTYPES
  private var log10PError = VariantContext.NO_LOG10_PERROR
  private var filters : Set[String] = null
  private var attributesMap : immutable.Map[String, Any] = null
  private var attributesCanBeModified = false

  /** enum of what must be validated */
  private var validationMode  = ValidationMode.ValueSet.empty




  /**
   * Create an empty VariantContextBuilder where all values adopt their default values, but the bare min.
   * of info (source, chr, start, stop, and alleles) have been provided to start.
   */
  @Requires(Array[String]("source != null", "contig != null", "start >= 0", "stop >= 0", "alleles != null && !alleles.isEmpty()"))
  def this( source : String,  contig : String,  start : Long,  stop : Long,  alleleContext : AlleleContext) {
    this()
    this.source = source
    this.contig = contig
    this.start = start
    this.stop = stop
    this.alleleContext = alleleContext
    this.attributesMap = immutable.HashMap[String, Any]()
    this.validationMode = ValidationMode.ValueSet.empty + ValidationMode.ALLELES
  }

  /**
   * Returns a new builder based on parent -- the new VC will have all fields initialized
   * to their corresponding values in parent.  This is the best way to create a derived VariantContext
   *
   * @param parent  Cannot be null
   */
  def this( parent : VariantContext ) {
    this()
    if ( parent == null ) throw new IllegalArgumentException("BUG: VariantContextBuilder parent argument cannot be null in VariantContextBuilder");
    this.alleleContext = parent.alleleContext
    this.attributesMap = parent.getAttributes
    this.attributesCanBeModified = false
    this.contig = parent.contig
    this.filters = parent.getFiltersMaybeNull
    this.genotypeContext = parent.genotypeContext
    this.ID = parent.getID
    this.log10PError = parent.getLog10PError
    this.source = parent.getSource
    this.start = parent.getStart()
    this.stop = parent.getEnd()
    this.fullyDecoded = parent.isFullyDecoded()
  }

  def this ( parent : VariantContextBuilder) {
    this()
    if ( parent == null ) throw new IllegalArgumentException("BUG: VariantContext parent argument cannot be null in VariantContextBuilder");
    this.alleleContext = parent.alleleContext
    this.attributesCanBeModified = false
    this.contig = parent.contig
    this.genotypeContext = parent.genotypeContext
    this.ID = parent.ID
    this.log10PError = parent.log10PError
    this.source = parent.source
    this.start = parent.start
    this.stop = parent.stop
    this.fullyDecoded = parent.fullyDecoded

    this.attributes(parent.attributesMap)
    this.filters(parent.filters)
  }

  def copy() : VariantContextBuilder = { new  VariantContextBuilder(this); }

  /**
   * Tells this builder to use this collection of alleles for the resulting VariantContext
   *
   * @param alleleContext
   * @return this builder
   */
  @Requires(Array[String]("alleles != null", "!alleles.isEmpty()"))
  def alleles(  alleleContext : AlleleContext) : VariantContextBuilder = {
    this.alleleContext = alleleContext
    validationMode += ValidationMode.ALLELES
    this
  }

  def alleles(  alleleStrings : Array[String]) : VariantContextBuilder =  {

    val alleleArray = new Array[Allele](alleleStrings.size)
    var counter = 0
    while(counter < alleleStrings.size)
    {
      alleleArray(counter) = Allele(alleleStrings(counter), counter == 0)
      counter += 1
    }
    val alleleContext = new AlleleContext(alleleArray)

    alleles(alleleContext);
  }

  def alleles(allelesArray : Array[Allele]): VariantContextBuilder ={


    alleles(new AlleleContext(allelesArray));


  }




  def alleles( alleleStrings : String*) : VariantContextBuilder = {
    val alleleStringArray = new Array[String](alleleStrings.size)
    var counter = 0
    while(counter < alleleStrings.size)
    {
      alleleStringArray(counter) = alleleStrings(counter)
      counter += 1
    }
    alleles(alleleStringArray)
  }

  def getAlleles() : AlleleContext = {  alleleContext  }

  /**
   * Tells this builder to use this map of attributes alleles for the resulting VariantContext
   *
   * Attributes can be null -> meaning there are no attributes.  After
   * calling this routine the builder assumes it can modify the attributes
   * object here, if subsequent calls are made to set attribute values
   * @param attributesMap
   */
  def attributes(  attributesMap : immutable.Map[String, Any]) : VariantContextBuilder=  {
    if (attributesMap != null) {
      this.attributesMap = attributesMap
    }
    else {
      this.attributesMap = immutable.HashMap[String, Any]()
    }

    this.attributesCanBeModified = true
    this
  }

  //passes on a mutable map from java as an immutable map
  def attributesFromJava(  attributesMap : Map[String, Any]) : VariantContextBuilder =
  {
    attributes(attributesMap.toMap)
  }

  /**
   * Puts the key -> value mapping into this builder's attributes
   *
   * @param key
   * @param value
   * @return
   */
  @Requires(Array[String]("key != null"))
  @Ensures(Array[String]("this.attributes.size() == old(this.attributes.size()) || this.attributes.size() == old(this.attributes.size()+1)"))
  def attribute(  key : String,  value : Any) : VariantContextBuilder = {
    attributesMap += (key -> value)
    this;
  }

  /**
   * Removes key if present in the attributes
   *
   * @param key  key to remove
   * @return
   */
  @Requires(Array[String]("key != null"))
  @Ensures(Array[String]("this.attributes.size() == old(this.attributes.size()) || this.attributes.size() == old(this.attributes.size()-1)"))
  def rmAttribute( key : String) : VariantContextBuilder = {
    attributesMap -= (key)
    this;
  }

  /**
   * Removes list of keys if present in the attributes
   *
   * @param keys  list of keys to remove
   * @return
   */
  @Requires(Array[String]("keys != null"))
  @Ensures(Array[String]("this.attributes.size() <= old(this.attributes.size())"))
  def rmAttributes( keys : Array[String]) : VariantContextBuilder = {
    attributesMap --= keys
    this
  }


  /**
   * This builder's filters are set to this value
   *
   * filters can be null -> meaning there are no filters
   * @param filters
   */
  def filters( filters : Set[String]) : VariantContextBuilder = {
    this.filters = filters
    this
  }

  /**
   * {@link #filters}
   *
   * @param filtersSeq
   * @return
   */
  def filters( filtersSeq : String*) : VariantContextBuilder =  {
    val linkedHashSet = new mutable.LinkedHashSet[String]()
    linkedHashSet ++= filtersSeq

    filters(linkedHashSet);
    this;
  }

  @Requires(Array[String]("filter != null", "!filter.equals(\"PASS\")"))
  def filter( filter : String) : VariantContextBuilder  = {
    if ( this.filters == null ){ this.filters = new mutable.LinkedHashSet[String]()}
    this.filters += (filter)
    this
  }

  /**
   * Tells this builder that the resulting VariantContext should have PASS filters
   *
   * @return
   */
  def passFilters() : VariantContextBuilder = {
    filters(VariantContext.PASSES_FILTERS);
  }

  /**
   * Tells this builder that the resulting VariantContext be unfiltered
   *
   * @return
   */
  def unfiltered() : VariantContextBuilder = {
    this.filters = null
    this
  }

  /**
   * Tells this builder that the resulting VariantContext should use this genotypes GenotypeContext
   *
   * Note that genotypes can be null -> meaning there are no genotypes
   *
   * @param genotypes
   */
  def genotypeContext(  genotypes : GenotypesContext) : VariantContextBuilder = {
    this.genotypeContext = genotypes;
    if ( genotypes != null )
      validationMode += ValidationMode.GENOTYPES
    this;
  }

  def genotypesNoValidation(  genotypes : GenotypesContext) : VariantContextBuilder = {
    this.genotypeContext = genotypes;
    this;
  }

  /**
   * Tells this builder that the resulting VariantContext should use a GenotypeContext containing genotypes
   *
   * Note that genotypes can be null -> meaning there are no genotypes
   *
   * @param genotypesSeq
   */
  def genotypes( genotypesSeq : Seq[Genotype]) : VariantContextBuilder = {
    genotypeContext(GenotypesContext.copy(genotypesSeq));
  }



  /**
   * Tells this builder that the resulting VariantContext should not contain any GenotypeContext
   */
  def noGenotypes() : VariantContextBuilder = {
    this.genotypeContext = null;
    this;
  }

  /**
   * Tells us that the resulting VariantContext should have ID
   * @param ID
   * @return
   */
  @Requires(Array[String]("ID != null"))
  def id( ID : String)  : VariantContextBuilder = {
    this.ID = ID;
    this;
  }

  /**
   * Tells us that the resulting VariantContext should not have an ID
   * @return
   */
  def noID() : VariantContextBuilder =  {
    id(VCFConstants.EMPTY_ID_FIELD);
  }

  /**
   * Tells us that the resulting VariantContext should have log10PError
   * @param log10PError
   * @return
   */
  @Requires(Array[String]("log10PError <= 0 || log10PError == VariantContext.NO_LOG10_PERROR"))
  def log10PError( log10PError : Double) : VariantContextBuilder ={
    this.log10PError = log10PError;
    this;
  }

  /**
   * Tells us that the resulting VariantContext should have source field set to source
   * @param source
   * @return
   */
  @Requires(Array[String]("source != null"))
  def source( source : String) : VariantContextBuilder = {
    this.source = source;
    this;
  }

  /**
   * Tells us that the resulting VariantContext should have the specified location
   * @param contig
   * @param start
   * @param stop
   * @return
   */
  @Requires(Array[String]("contig != null", "start >= 0", "stop >= 0"))
  def loc(  contig : String , start : Long, stop : Long )  : VariantContextBuilder = {
    this.contig = contig;
    this.start = start;
    this.stop = stop;
    validationMode += ValidationMode.ALLELES
    this;
  }

  /**
   * Tells us that the resulting VariantContext should have the specified contig chr
   * @param contig
   * @return
   */
  @Requires(Array[String]("contig != null"))
  def chr( contig : String)  :VariantContextBuilder ={
    this.contig = contig;
    this;
  }

  /**
   * Tells us that the resulting VariantContext should have the specified contig start
   * @param start
   * @return
   */
  @Requires(Array[String]("start >= 0"))
  def start(  start : Long ) : VariantContextBuilder = {
    this.start = start;
    validationMode += ValidationMode.ALLELES
    this;
  }

  /**
   * Tells us that the resulting VariantContext should have the specified contig stop
   * @param stop
   * @return
   */
  @Requires(Array[String]("stop >= 0"))
  def  stop( stop : Long) : VariantContextBuilder = {
    this.stop = stop;
    this;
  }

  /**
   * @see #computeEndFromAlleles(java.utils.List, int, int) with endForSymbolicAlleles == -1
   */
  def computeEndFromAlleles(  alleles : Array[Allele] , start : Int) : VariantContextBuilder = {
    computeEndFromAlleles(alleles, start, -1);
  }

  /**
   * Compute the end position for this VariantContext from the alleles themselves
   *
   * assigns this builder the stop position computed.
   *
   * @param alleles the list of alleles to consider.  The reference allele must be the first one
   * @param start the known start position of this event
   * @param endForSymbolicAlleles the end position to use if any of the alleles is symbolic.  Can be -1
   *                              if no is expected but will throw an error if one is found
   * @return this builder
   */
  @Requires(Array[String]("! alleles.isEmpty()", "start > 0", "endForSymbolicAlleles == -1 || endForSymbolicAlleles > 0" ))
  def computeEndFromAlleles(  alleles : Array[Allele] ,  start :  Int,  endForSymbolicAlleles : Int) : VariantContextBuilder = {
    stop(VariantContextUtils.computeEndFromAlleles(alleles, start, endForSymbolicAlleles))
    this;
  }

  /**
   * @return true if this builder contains fully decoded data
   *
   * See VariantContext for more information
   */
  def  isFullyDecoded() : Boolean  = { fullyDecoded }

  /**
   * Sets this builder's fully decoded state to true.
   *
   * A fully decoded builder indicates that all fields are represented by their
   * proper java objects (e.g., Integer(10) not "10").
   *
   * See VariantContext for more information
   *
   * @param isFullyDecoded
   */
  def fullyDecoded( isFullyDecoded : Boolean )  : VariantContextBuilder = {
    this.fullyDecoded = isFullyDecoded;
    this;
  }

  /**
   * Takes all of the builder data provided up to this point, and instantiates
   * a freshly allocated VariantContext with all of the builder data.  This
   * VariantContext is validated as appropriate and if not failing QC (and
   * throwing an exception) is returned.
   *
   * Note that this function can be called multiple times to create multiple
   * VariantContexts from the same builder.
   */
  def  make() : VariantContext = {
    new VariantContext(source, ID, contig, start, stop, alleleContext,
      genotypeContext, log10PError, filters, attributesMap,
      fullyDecoded, validationMode);
  }


}

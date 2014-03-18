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

/**
 * Class VariantContext
 *
 * == High-level overview ==
 *
 * The VariantContext object is a single general class system for representing genetic variation data composed of:
 *
 * * Allele: representing single genetic haplotypes (A, T, ATC, -) (note that null alleles are used here for illustration; see the Allele class for how to represent indels)
 * * Genotype: an assignment of alleles for each chromosome of a single named sample at a particular locus
 * * VariantContext: an abstract class holding all segregating alleles at a locus as well as genotypes
 *    for multiple individuals containing alleles at that locus
 *
 * The class system works by defining segregating alleles, creating a variant context representing the segregating
 * information at a locus, and potentially creating and associating genotypes with individuals in the context.
 *
 * All of the classes are highly validating -- call validate() if you modify them -- so you can rely on the
 * self-consistency of the data once you have a VariantContext in hand.  The system has a rich set of assessor
 * and manipulator routines, as well as more complex static support routines in VariantContextUtils.
 *
 * The VariantContext (and Genotype) objects are attributed (supporting addition of arbitrary key/value pairs) and
 * filtered (can represent a variation that is viewed as suspect).
 *
 * VariantContexts are dynamically typed, so whether a VariantContext is a SNP, Indel, or NoVariant depends
 * on the properties of the alleles in the context.  See the detailed documentation on the Type parameter below.
 *
 * It's also easy to create subcontexts based on selected genotypes.
 *
 * == Working with Variant Contexts ==
 * By default, VariantContexts are immutable.  In order to access (in the rare circumstances where you need them)
 * setter routines, you need to create MutableVariantContexts and MutableGenotypes.
 *
 * === Some example data ===
 *
 * Allele A, Aref, T, Tref;
 * Allele del, delRef, ATC, ATCref;
 *
 * A [ref] / T at 10
 * GenomeLoc snpLoc = GenomeLocParser.createGenomeLoc("chr1", 10, 10);
 *
 * A / ATC [ref] from 20-23
 * GenomeLoc delLoc = GenomeLocParser.createGenomeLoc("chr1", 20, 22);
 *
 *  // A [ref] / ATC immediately after 20
 * GenomeLoc insLoc = GenomeLocParser.createGenomeLoc("chr1", 20, 20);
 *
 * === Alleles ===
 *
 * See the documentation in the Allele class itself
 *
 * What are they?
 *
 * Alleles can be either reference or non-reference
 *
 * Examples of alleles used here:
 *
 *   A = new Allele("A");
 *   Aref = new Allele("A", true);
 *   T = new Allele("T");
 *   ATC = new Allele("ATC");
 *
 * === Creating variant contexts ===
 *
 * ==== By hand ====
 *
 * Here's an example of a A/T polymorphism with the A being reference:
 *
 * <pre>
 * VariantContext vc = new VariantContext(name, snpLoc, Arrays.asList(Aref, T));
 * </pre>
 *
 * If you want to create a non-variant site, just put in a single reference allele
 *
 * <pre>
 * VariantContext vc = new VariantContext(name, snpLoc, Arrays.asList(Aref));
 * </pre>
 *
 * A deletion is just as easy:
 *
 * <pre>
 * VariantContext vc = new VariantContext(name, delLoc, Arrays.asList(ATCref, del));
 * </pre>
 *
 * The only thing that distinguishes between an insertion and deletion is which is the reference allele.
 * An insertion has a reference allele that is smaller than the non-reference allele, and vice versa for deletions.
 *
 * <pre>
 * VariantContext vc = new VariantContext("name", insLoc, Arrays.asList(delRef, ATC));
 * </pre>
 *
 * ==== Converting rods and other data structures to VCs ====
 *
 * You can convert many common types into VariantContexts using the general function:
 *
 * <pre>
 * VariantContextAdaptors.convertToVariantContext(name, myObject)
 * </pre>
 *
 * dbSNP and VCFs, for example, can be passed in as myObject and a VariantContext corresponding to that
 * object will be returned.  A null return type indicates that the type isn't yet supported.  This is the best
 * and easiest way to create contexts using RODs.
 *
 *
 * === Working with genotypes ===
 *
 * <pre>
 * List<Allele> alleles = Arrays.asList(Aref, T);
 * Genotype g1 = new Genotype(Arrays.asList(Aref, Aref), "g1", 10);
 * Genotype g2 = new Genotype(Arrays.asList(Aref, T), "g2", 10);
 * Genotype g3 = new Genotype(Arrays.asList(T, T), "g3", 10);
 * VariantContext vc = new VariantContext(snpLoc, alleles, Arrays.asList(g1, g2, g3));
 * </pre>
 *
 * At this point we have 3 genotypes in our context, g1-g3.
 *
 * You can assess a good deal of information about the genotypes through the VariantContext:
 *
 * <pre>
 * vc.hasGenotypes()
 * vc.isMonomorphicInSamples()
 * vc.isPolymorphicInSamples()
 * vc.getSamples().size()
 *
 * vc.getGenotypes()
 * vc.getGenotypes().get("g1")
 * vc.hasGenotype("g1")
 *
 * vc.getCalledChrCount()
 * vc.getCalledChrCount(Aref)
 * vc.getCalledChrCount(T)
 * </pre>
 *
 * === NO_CALL alleles ===
 *
 * The system allows one to create Genotypes carrying special NO_CALL alleles that aren't present in the
 * set of context alleles and that represent undetermined alleles in a genotype:
 *
 * Genotype g4 = new Genotype(Arrays.asList(Allele.NO_CALL, Allele.NO_CALL), "NO_DATA_FOR_SAMPLE", 10);
 *
 *
 * === subcontexts ===
 * It's also very easy get subcontext based only the data in a subset of the genotypes:
 *
 * <pre>
 * VariantContext vc12 = vc.subContextFromGenotypes(Arrays.asList(g1,g2));
 * VariantContext vc1 = vc.subContextFromGenotypes(Arrays.asList(g1));
 * </pre>
 *
 * <s3>
 *     Fully decoding.  Currently VariantContexts support some fields, particularly those
 *     stored as generic attributes, to be of any type.  For example, a field AB might
 *     be naturally a floating point number, 0.51, but when it's read into a VC its
 *     not decoded into the Java presentation but left as a string "0.51".  A fully
 *     decoded VariantContext is one where all values have been converted to their
 *     corresponding Java object types, based on the types declared in a VCFHeader.
 *
 *     The fullyDecode() takes a header object and creates a new fully decoded VariantContext
 *     where all fields are converted to their true java representation.  The VCBuilder
 *     can be told that all fields are fully decoded, in which case no work is done when
 *     asking for a fully decoded version of the VC.
 * </s3>
 *
 * @author depristo
 */

/**
 * Created by Wim Spee on 2/20/14.
 */




import org.broad.tribble.util.ParsingUtils;
import org.broadinstitute.variant.utils.GeneralUtils;
import org.broadinstitute.variant.vcf.VCFConstants
import scala.collection._
import scala.collection.JavaConversions._
import org.broadinstitute.variant.vcf._
import org.broad.tribble.{TribbleException, Feature}


// ---------------------------------------------------------------------------------------------------------
//
// constructors: see VariantContextBuilder
//
// ---------------------------------------------------------------------------------------------------------

/**
 * the actual constructor.  Private access only
 *
 * @param source          source
 * @param contig          the contig
 * @param start           the start base (one based)
 * @param stop            the stop reference base (one based)
 * @param alleleContext         alleles
 * @param _genotypeContext       genotypes map
 * @param log10PError  qual
 * @param filters         filters: use null for unfiltered and empty set for passes filters
 * @param attributes      attributes
 * @param validationToPerform     set of validation steps to take
 */
class VariantContext(
                      val source : String,
                      val _ID : String,
                      val contig : String,
                      val start : Long,
                      val stop : Long,
                      val alleleContext : AlleleContext,
                      private val _genotypeContext : GenotypesContext,
                      val log10PError : Double,
                      val filters : Set[String],
                      val attributes : immutable.Map[String, Any],
                      val fullyDecoded : Boolean,
                      val validationToPerform : ValidationMode.ValueSet
                      ) extends Feature{  // to enable tribble integration

  protected val commonInfo : CommonInfo = new CommonInfo(source, log10PError, filters, attributes)

  if ( contig == null ) { throw new IllegalArgumentException("Contig cannot be null"); }

  // intern for efficiency.  equals calls will generate NPE if ID is inappropriately passed in as null
  if ( _ID == null || _ID.equals("") ){ throw new IllegalArgumentException("ID field cannot be the null or the empty string");}
  val ID = if( _ID.equals(VCFConstants.EMPTY_ID_FIELD)) {VCFConstants.EMPTY_ID_FIELD}else{ _ID}

  /** A mapping from sampleName -> genotype objects for all genotypes associated with this context */
  val genotypeContext : GenotypesContext = if ( _genotypeContext == null || _genotypeContext == VariantContext.NO_GENOTYPES ){  VariantContext.NO_GENOTYPES;  }
  else { _genotypeContext.immutable();   }

  private var mNotYetComputed = true
  private var monomorphic : Boolean = false;            // cached monomorphic value: null -> not yet computed, False, True

  validate(validationToPerform)


  /**
   * Copy constructor
   *
   * @param other the VariantContext to copy
   */
  protected def this( other : VariantContext) {
    this(other.getSource, other.getID, other.getChr, other.getStart, other.getEnd,
      other.alleleContext, other.genotypeContext, other.getLog10PError,
      other.getFiltersMaybeNull,
      other.getAttributes,
      other.fullyDecoded, VariantContext.NO_VALIDATION);
  }



  // ---------------------------------------------------------------------------
  //
  // Monomorphic methods
  //
  // ---------------------------------------------------------------------------

  /**
   * First check if we are not a an variant according to AlleleContext. If we are a variant check if we are monomorphic in the samples.
   * Genotype-specific functions -- are the genotypes monomorphic w.r.t. to the alleles segregating at this
   * site?  That is, is the number of alternate alleles among all fo the genotype == 0?
   *
   * @return true if it's monomorphic
   */
  def isMonomorphic() : Boolean = {
    if ( mNotYetComputed){
      monomorphic = ! alleleContext.isVariant() || genotypeContext.isMonomorphicInSamples(alleleContext.getReference())
      mNotYetComputed = false
    }
    monomorphic;
  }

  /**
   * Genotype-specific functions -- are the genotypes polymorphic w.r.t. to the alleles segregating at this
   * site?  That is, is the number of alternate alleles among all fo the genotype > 0?
   *
   * @return true if it's polymorphic
   */
  def isPolymorphic(): Boolean = { ! isMonomorphic();  }

  /**
   * @return true if the context has associated genotypes
   */
  def hasGenotypes: Boolean = { !genotypeContext.isEmpty }

  /*
  * Determine which genotype fields are in use in the genotypes in VC
  * @return an ordered list of genotype fields in use in VC.  If vc has genotypes this will always include GT first
  */
  def calcVCFGenotypeKeys( header : VCFHeader) : Array[String] = {

    val keys = scala.collection.mutable.HashSet[String]();



    var sawGoodGT = false;
    var sawGoodQual = false;
    var sawGenotypeFilter = false;
    var sawDP = false;
    var sawAD = false;
    var sawPL = false;

    var counter = 0
    val limit = genotypeContext.size

    while(counter < limit)
    {
      val g = genotypeContext(counter)
      keys ++= g.getExtendedAttributes().keySet
      if ( g.isAvailable() ) sawGoodGT = true;
      if ( g.hasGQ() ) sawGoodQual = true;
      if ( g.hasDP() ) sawDP = true;
      if ( g.hasAD() ) sawAD = true;
      if ( g.hasPL() ) sawPL = true;
      if (g.isFiltered()) sawGenotypeFilter = true;
      counter +=1
    }

    if ( sawGoodQual ) keys.add(VCFConstants.GENOTYPE_QUALITY_KEY);
    if ( sawDP ) keys.add(VCFConstants.DEPTH_KEY);
    if ( sawAD ) keys.add(VCFConstants.GENOTYPE_ALLELE_DEPTHS);
    if ( sawPL ) keys.add(VCFConstants.GENOTYPE_PL_KEY);
    if ( sawGenotypeFilter ) keys.add(VCFConstants.GENOTYPE_FILTER_KEY);

    var sortedList = keys.toArray.sorted
    //var sortedList = ParsingUtils.sortList(new ArrayList[String](keys));

    // make sure the GT is first
    if (sawGoodGT) { sortedList = VCFConstants.GENOTYPE_KEY +: sortedList }

    if (sortedList.isEmpty && header.hasGenotypingData) {
      // this needs to be done in case all samples are no-calls
      Array[String](VCFConstants.GENOTYPE_KEY);
    } else {
      sortedList;
    }
  }


  // ---------------------------------------------------------------------------------------------------------
  //
  // Selectors
  //
  // ---------------------------------------------------------------------------------------------------------

  /**
   * This method subsets down to a set of samples.
   *
   * At the same time returns the alleles to just those in use by the samples,
   * if rederiveAllelesFromGenotypes is true, otherwise the full set of alleles
   * in this VC is returned as the set of alleles in the subContext, even if
   * some of those alleles aren't in the samples
   *
   * WARNING: BE CAREFUL WITH rederiveAllelesFromGenotypes UNLESS YOU KNOW WHAT YOU ARE DOING?
   *
   * @param sampleNames    the sample names
   * @param rederiveAllelesFromGenotypes if true, returns the alleles to just those in use by the samples, true should be default
   * @return new VariantContext subsetting to just the given samples
   */
  def subContextFromSamples( sampleNames : scala.collection.Set[String],  rederiveAllelesFromGenotypes : Boolean ) : VariantContext = {
    if ( sampleNames.containsAll(genotypeContext.getSampleNames()) && ! rederiveAllelesFromGenotypes ) {
      this; // fast path when you don't have any work to do
    } else {
      val builder = new VariantContextBuilder(this);
      val newGenotypes = genotypeContext.subsetToSamples(sampleNames);

      if ( rederiveAllelesFromGenotypes ) {
        val allelesFromGenotypes = allelesOfGenotypes(newGenotypes);

        // ensure original order of genotypes
        val rederivedAlleles = alleleContext.getAlleles().filter(allelesFromGenotypes.contains(_))
        builder.alleles(rederivedAlleles);
      }
      else {
        builder.alleles(alleleContext.getAlleles());
      }

      builder.genotypes(newGenotypes).make();
    }
  }

  /**
   * @see #subContextFromSamples(java.utils.Set, boolean) with rederiveAllelesFromGenotypes = true
   *
   * @param sampleNames
   * @return
   */
  def subContextFromSamples( sampleNames : Set[String])  : VariantContext = { subContextFromSamples(sampleNames, true);             }
  def subContextFromSample( sampleName : String)  : VariantContext =        { subContextFromSamples(immutable.HashSet(sampleName)); }

  /**
   * helper routine for subcontext
   * @param genotypes genotypes
   * @return allele set
   */
  private def allelesOfGenotypes( genotypes : Seq[Genotype]) :Set[Allele] = {

    val calledAllelesArray = genotypes.flatMap(_.getAlleles()).filter(_.isCalled).distinct

    if(calledAllelesArray.count(_.isReference) == 0)
    {
      (calledAllelesArray :+ alleleContext.getReference()).toSet
    }
    else
    {
      calledAllelesArray.toSet
    }
  }


  // ---------------------------------------------------------------------------------------------------------
  //
  // Generic accessors
  //
  // ---------------------------------------------------------------------------------------------------------

  def getID: String     = { ID;                                     }
  def hasID : Boolean   = { getID != VCFConstants.EMPTY_ID_FIELD;   }
  def emptyID :Boolean  = { ! hasID;                                }




  // ---------------------------------------------------------------------------------------------------------
  //
  // get routines to access context info fields
  //
  // ---------------------------------------------------------------------------------------------------------
  def getSource             = {  commonInfo.getName }
  def getFiltersMaybeNull   = {  commonInfo.getFiltersMaybeNull }
  def getFilters            = {  commonInfo.getFilters }
  def isFiltered            = {  commonInfo.isFiltered }
  def isNotFiltered         = {  commonInfo.isNotFiltered }
  def filtersWereApplied    = {  commonInfo.filtersWereApplied }
  def hasLog10PError        = {  commonInfo.hasLog10PError }
  def getLog10PError        = {  commonInfo.getLog10PError }
  def getPhredScaledQual    = {  commonInfo.getPhredScaledQual }

  def getAttributes :  immutable.Map[String, Any]   = {  commonInfo.getAttributes(); }
  def hasAttribute( key : String) : Boolean = {  commonInfo.hasAttribute(key); }
  def getAttribute( key : String) : Any     = {  commonInfo.getAttribute(key); }

  def getAttribute( key : String, defaultValue : Any) : Any = {
    return commonInfo.getAttribute(key, defaultValue);
  }

  def getAttributeAsString( key : String,  defaultValue :String)  = { commonInfo.getAttributeAsString(key, defaultValue); }
  def getAttributeAsInt( key :String,  defaultValue :Int)         = { commonInfo.getAttributeAsInt(key, defaultValue); }
  def getAttributeAsDouble( key : String,   defaultValue :Double) = { commonInfo.getAttributeAsDouble(key, defaultValue); }
  def getAttributeAsBoolean( key :String, defaultValue : Boolean) = { commonInfo.getAttributeAsBoolean(key, defaultValue); }

  def getCommonInfo() = { commonInfo;   }



  /**
   * Returns the maximum ploidy of all samples in this VC, or default if there are no genotypes
   *
   * This function is caching, so it's only expensive on the first call
   *
   * @param defaultPloidy the default ploidy, if all samples are no-called
   * @return default, or the max ploidy
   */
  def getMaxPloidy(defaultPloidy : Int) = { genotypeContext.getMaxPloidy(defaultPloidy);  }


  // ---------------------------------------------------------------------------------------------------------
  //
  // validation: extra-strict validation routines for paranoid users
  //
  // ---------------------------------------------------------------------------------------------------------

  /**
   * Run all extra-strict validation tests on a Variant Context object
   *
   * @param reportedReference   the reported reference allele
   * @param observedReference   the actual reference allele
   * @param rsIDs               the true dbSNP IDs
   */
  def extraStrictValidation(  reportedReference : Allele,   observedReference : Allele,   rsIDs : Set[String]) {
    // validate the reference
    validateReferenceBases(reportedReference, observedReference);

    // validate the RS IDs
    validateRSIDs(rsIDs);

    // validate the altenate alleles
    validateAlternateAlleles();

    // validate the AN and AC fields
    validateChromosomeCounts();

    // TODO: implement me
    //checkReferenceTrack();
  }

  def validateReferenceBases(  reportedReference : Allele,  observedReference : Allele) {
    if ( reportedReference != null && !reportedReference.basesMatch(observedReference) ) {
      val position =  getChr()+":" + getStart().toString; val observedRefBaseString = observedReference.getBaseString(); val reportedRefBaseString = reportedReference.getBaseString()
      throw new TribbleException.InternalCodecException(s"the REF allele is incorrect for the record at position $position , fasta says $observedRefBaseString vs. VCF says $reportedRefBaseString");
    }
  }

  def validateRSIDs( rsIDs : Set[String]) {
    if ( rsIDs != null && hasID ) {
      val splitIds =  getID.split(VCFConstants.ID_FIELD_SEPARATOR)

      var counter = 0
      val limit = splitIds.size

      while(counter < limit)
      {
        val id = splitIds(counter)
        if ( id.startsWith("rs") && !rsIDs.contains(id) ){
          val position =  getChr()+":" + getStart().toString;
          throw new TribbleException.InternalCodecException(s"the rsID $id for the record at position $position is not in dbSNP")}
        counter +=1
      }
    }
  }

  /**
   * Validate that the AlleleContext and the GenotypeContext contain the same alternate alleles
   *
   */
  def validateAlternateAlleles() {
    if ( !hasGenotypes ){}
    else
    {
      val alleleContextAlleles = alleleContext.getAlleles();
      val genotypeContextAlleles = mutable.HashSet[Allele]();
      genotypeContextAlleles.add(alleleContext.getReference());

      genotypeContext.getGenotypes.filter(_.isCalled()).foreach(genotypeContextAlleles ++= _.getAlleles())

      if(genotypeContextAlleles.contains(Allele.NO_CALL)){genotypeContextAlleles -= Allele.NO_CALL}

      if ( alleleContextAlleles.size != genotypeContextAlleles.size ){
        val position =  getChr()+":" + getStart().toString
        throw new TribbleException.InternalCodecException(s"one or more of the ALT allele(s) for the record at position $position are not observed at all in the sample genotypes")}

      // take the intersection and see if things change
      val interSect = genotypeContextAlleles.filter(alleleContextAlleles.contains(_));

      if ( interSect.size != alleleContextAlleles.size ){
        val position =  getChr()+":" + getStart().toString
        throw new TribbleException.InternalCodecException(s"one or more of the ALT allele(s) for the record at position $position are not observed at all in the sample genotypes")}
    }
  }

  /**
   * Validate that the AN (alleleNumber) is equal as reported the variantcontext attributes and as observed in the GenotypeContext
   *
   */
  def validateAlleleNumber()
  {
    if ( hasAttribute(VCFConstants.ALLELE_NUMBER_KEY) ) {
      val position =  getChr()+":" + getStart().toString;
      val reportedAN = getAttributeAsInt(VCFConstants.ALLELE_NUMBER_KEY, 0)
      val genotypeContextAN = genotypeContext.getCalledChrCount();
      if ( reportedAN != genotypeContextAN ){
        throw new TribbleException.InternalCodecException(s"the Allele Number (AN) tag is incorrect for the record at position $position, $reportedAN vs. $genotypeContextAN")}
    }
  }

  def validateAlleleCount()
  {
    if ( hasAttribute(VCFConstants.ALLELE_COUNT_KEY) )
    {
      val position =  getChr()+":" + getStart().toString;

      val alternateAlleles = alleleContext.getAlternateAlleles()
      val genotypeContextACs = if(alternateAlleles.size > 0)
      {
        alternateAlleles.map(genotypeContext.getCalledChrCount(_))
      }
      else
      {
        Array[Int](0)
      }

      if( getAttribute(VCFConstants.ALLELE_COUNT_KEY).isInstanceOf[Array[Any]] ){ validateAlleleCountArray(genotypeContextACs) }
      else
      {
        if ( genotypeContextACs.size != 1 ){
          throw new TribbleException.InternalCodecException(s"the Allele Count (AC) tag doesn't have enough values for the record at position $position")}
        val reportedAC = getAttribute(VCFConstants.ALLELE_COUNT_KEY).asInstanceOf[String].toInt;
        val genotypeAC = genotypeContextACs(0)
        if ( reportedAC != genotypeAC ){
          throw new TribbleException.InternalCodecException(s"the Allele Count (AC) tag is incorrect for the record at position $position $reportedAC vs $genotypeAC")}
      }
    }
  }

  /**
   * Validate that the  AC(alleleCount) are equal as reported in variantcontext attributes and as observed in  GenotypeContext
   *
   */
  def validateAlleleCountArray(genotypeContextACs : Array[Int]) {

    val position =  getChr()+":" + getStart().toString;
    val reportedACs = getAttribute(VCFConstants.ALLELE_COUNT_KEY).asInstanceOf[Array[Any]];

    if (genotypeContextACs.size != reportedACs.size) {
      val reportedACCount = reportedACs.size;
      val genotypeACCount = genotypeContextACs.size
      throw new TribbleException.InternalCodecException(s"the Allele Count (AC) tag doesn't have the correct number of values for the record at position $position reported=$reportedACCount observed=$genotypeACCount")
    }
    else
    {
      val reportedACsasIntegers = reportedACs.map(_.asInstanceOf[String].toInt)
      var counter = 0
      val limit = genotypeContextACs.size
      while (counter < limit) {
        if (reportedACsasIntegers(counter) != genotypeContextACs(counter)) {
          val reportedACCount = reportedACsasIntegers(counter);
          val genotypeACCount = genotypeContextACs(counter)
          throw new TribbleException.InternalCodecException(s"the Allele Count (AC) tag is incorrect for the record at position $position, $reportedACCount vs. $genotypeACCount at allele index $counter");
        }
        counter += 1
      }
    }
  }


  /**
   * Validate that the AN (alleleNumber) and AC(alleleCount) are equal in the variantcontext attributes and the GenotypeContext
   *
   */
  def validateChromosomeCounts() {
    if ( hasGenotypes )
    {
      validateAlleleNumber()  // AN
      validateAlleleCount()   // AC
    }
  }

  // ---------------------------------------------------------------------------------------------------------
  //
  // validation: the normal validation routines are called automatically upon creation of the VC
  //
  // ---------------------------------------------------------------------------------------------------------

  def validate( validationToPerform : ValidationMode.ValueSet) : Boolean = {
    validateStop();

    if(validationToPerform.contains(ValidationMode.ALLELES)){ alleleContext.validateAlleles() }
    if(validationToPerform.contains(ValidationMode.GENOTYPES)){ validateGenotypes() }

    true;
  }

  /**
   * Check that getEnd() == END from the info field, if it's present
   */
  def validateStop() {
    if ( hasAttribute(VCFConstants.END_KEY) )
    {
      val end = getAttributeAsInt(VCFConstants.END_KEY, -1);
      //assert end != -1;
      if ( end != getEnd() )
      {
        val message = "Badly formed variant context at location " + getChr() + ":" +
          getStart() + "; getEnd() was " + getEnd() +
          " but this VariantContext contains an END key with value " + end;
        if ( GeneralUtils.DEBUG_MODE_ENABLED && VariantContext.WARN_ABOUT_BAD_END ) {
          System.err.println(message);
        }
        else {
          throw new TribbleException(message);
        }
      }
    } else {
      val length : Long= (stop - start) + 1;
      if ( !alleleContext.hasSymbolicAlleles && length != alleleContext.getReference().length() ) {
        throw new IllegalStateException("BUG: GenomeLoc " + contig + ":" + start + "-" + stop + " has a size == " + length + " but the variation reference allele has length " + alleleContext.getReference().length() + " this = " + this);
      }
    }
  }

  def validateGenotypes() {
    if ( this.genotypeContext == null ){ throw new IllegalStateException("Genotypes is null");}

    genotypeContext.filter(_.isAvailable()).flatMap(_.getAlleles()).distinct.foreach
    {
      gAllele =>     if(!alleleContext.hasAllele(gAllele) && gAllele.isCalled  ){throw new IllegalStateException("Allele in genotype " + gAllele + " not in the alleleContext " + alleleContext.getAlleles())}
    }
  }


  // ---------------------------------------------------------------------------------------------------------
  //
  // utility routines
  //
  // ---------------------------------------------------------------------------------------------------------


  override def  toString() : String = {
    // Note: passing genotypes to String.format() will implicitly decode the genotypes
    // This may not be desirable, so don't decode by default

    if(genotypeContext.isLazyWithData) {toStringUnparsedGenotypes()} else{ toStringDecodeGenotypes();}
  }

  def toStringDecodeGenotypes() : String = {
    val source = getSource
    val location = contig + ":" + (if(start - stop == 0){ start} else{start + "-" + stop})
    val phred = getPhredScaledQual
    val qual = if(hasLog10PError){ val phred = getPhredScaledQual; f"$phred%.2f"}else{ "."}
    val variantType = alleleContext.getType().toString
    val alleles = alleleContext.getAlleles().sorted.mkString(",")
    val attributes = ParsingUtils.sortedString(this.getAttributes)
    val genotypes = genotypeContext.getGenotypes.mkString(",")

    s"[VC $source @ $location Q$qual of type=$variantType alleles=$alleles attr=$attributes GT=$genotypes"
  }

  def toStringUnparsedGenotypes() : String =  {
    val source = getSource
    val location = contig + ":" + (if(start - stop == 0){ start} else{start + "-" + stop})
    val qual = if(hasLog10PError){ val phred = getPhredScaledQual; f"$phred%.2f"}else{ "."}
    val variantType = alleleContext.getType().toString
    val alleles = alleleContext.getAlleles().sorted.mkString(",")
    val attributes = ParsingUtils.sortedString(this.getAttributes)
    val genotypes = genotypeContext.getUnparsedGenotypeData

    s"[VC $source @ $location Q$qual of type=$variantType alleles=$alleles attr=$attributes GT=$genotypes"

  }

  def toStringWithoutGenotypes() : String = {
    val source = getSource
    val location = contig + ":" + (if(start - stop == 0){ start} else{start + "-" + stop})
    val qual = if(hasLog10PError){ val phred = getPhredScaledQual; f"$phred%.2f"}else{ "."}
    val variantType = alleleContext.getType().toString
    val alleles = alleleContext.getAlleles().sorted.mkString(",")
    val attributes = ParsingUtils.sortedString(this.getAttributes)

    s"[VC $source @ $location Q$qual of type=$variantType alleles=$alleles attr=$attributes"
  }



  // ---------------------------------------------------------------------------------------------------------
  //
  // Fully decode
  //
  // ---------------------------------------------------------------------------------------------------------

  /**
   * Return a VC equivalent to this one but where all fields are fully decoded
   *
   * See VariantContext document about fully decoded
   *
   * @param header containing types about all fields in this VC
   * @return a fully decoded version of this VC
   */
  def fullyDecode(  header : VCFHeader, lenientDecoding : Boolean)  : VariantContext = {
    if ( isFullyDecoded() ){ this;}
    else {
      // TODO -- warning this is potentially very expensive as it creates copies over and over
      val builder = new VariantContextBuilder(this);
      fullyDecodeInfo(builder, header, lenientDecoding);
      fullyDecodeGenotypes(builder, header);
      builder.fullyDecoded(true);
      builder.make();
    }
  }

  /**
   * See VariantContext document about fully decoded
   * @return true if this is a fully decoded VC
   */
  def isFullyDecoded() = { fullyDecoded;  }

  private def fullyDecodeInfo( builder : VariantContextBuilder,  header : VCFHeader,  lenientDecoding : Boolean ) {
    builder.attributes(fullyDecodeAttributes(getAttributes, header, lenientDecoding));
  }

  private def fullyDecodeAttributes( attributes : immutable.Map[String, Any], header : VCFHeader, lenientDecoding : Boolean) : immutable.Map[String, Any] =
  {
    var newAttributes = immutable.HashMap[String, Any]()

    val attributesEntrySetIter =   attributes.entrySet().iterator()
    var counter = 0
    val limit = attributes.size

    while(counter < limit)
    {
      val attr = attributesEntrySetIter.next()
      val field = attr.getKey();

      if ( field.equals(VCFConstants.GENOTYPE_FILTER_KEY) ){} // gross, FT is part of the extended attributes
      else
      {
        val format = VariantContextUtils.getMetaDataForField(header, field);
        val decoded = decodeValue(field, attr.getValue(), format);

        if ( decoded != null && !lenientDecoding && format.getCountType() != VCFHeaderLineCount.UNBOUNDED  && format.getType() != VCFHeaderLineType.Flag )
        { // we expect exactly the right number of elements
        val obsSize = if(decoded.isInstanceOf[Array[Any]]){ decoded.asInstanceOf[Array[Any]].size} else{ 1 }
          val expSize = format.getCount(this);
          if ( obsSize != expSize ) {
            throw new TribbleException.InvalidHeader("Discordant field size detected for field " +
              field + " at " + getChr() + ":" + getStart() + ".  Field had " + obsSize + " values " +
              "but the header says this should have " + expSize + " values based on header record " +
              format);
          }
        }
        newAttributes += (field -> decoded)
      }
      counter +=1
    }
    newAttributes;
  }

  private def decodeValue( field : String,  value : Any,  format : VCFCompoundHeaderLine) : Any = {
    if ( value.isInstanceOf[String] )
    {
      if ( field.equals(VCFConstants.GENOTYPE_PL_KEY) ){ GenotypeLikelihoods.fromPLField(value.asInstanceOf[String]);}
      else
      {
        val string = value.asInstanceOf[String];
        if ( string.indexOf(",") != -1 ) {
          string.split(",").map(splitString => decodeOne(field, splitString, format));
        }
        else { decodeOne(field, string, format);  }
      }
    }
    else if ( value.isInstanceOf[Array[Any]] && (( value.asInstanceOf[Array[Any]](0)).isInstanceOf[String] ) )
    {
      val asList = value.asInstanceOf[Array[String]];
      asList.map(item => decodeOne(field, item, format))
    }
    else { value; }

    // allowMissingValuesComparedToHeader
  }

  private def decodeOne( field : String, string : String,  format : VCFCompoundHeaderLine) : Any = {
    try {
      if ( string.equals(VCFConstants.MISSING_VALUE_v4) ){ null }
      else
      {
        ( format.getType() ) match {
          case VCFHeaderLineType.Character => string;
          case VCFHeaderLineType.Flag =>
            val b = (string.toBoolean) || string.equals("1");
            if ( b == false ){
              throw new TribbleException("VariantContext FLAG fields " + field + " cannot contain false values"
                + " as seen at " + getChr() + ":" + getStart());}
            b;
          case VCFHeaderLineType.String   => string;
          case VCFHeaderLineType.Integer  => string.toInt
          case VCFHeaderLineType.Float    => string.toDouble;
          case _        => throw new TribbleException("Unexpected type for field" + field);
        }
      }
    } catch {
      case e : NumberFormatException =>  throw new TribbleException("Could not decode field " + field + " with value " + string + " of declared type " + format.getType());
    }
  }

  private def fullyDecodeGenotypes(  builder : VariantContextBuilder,  header : VCFHeader)   {
    val gc =  GenotypesContext.create();

    genotypeContext.getGenotypes.foreach(g => gc += (fullyDecodeGenotypes(g, header)))
    builder.genotypesNoValidation(gc);
  }

  private def fullyDecodeGenotypes(  g : Genotype,   header : VCFHeader) : Genotype = {
    val map = fullyDecodeAttributes(g.getExtendedAttributes(), header, true);
    new GenotypeBuilder(g).attributes(map).make();
  }

  // ---------------------------------------------------------------------------------------------------------
  //
  // tribble integration routines -- not for public consumption
  //
  // ---------------------------------------------------------------------------------------------------------
  override def getChr   = { contig       }
  override def getStart = { start.toInt  }
  override def getEnd   = { stop.toInt   }





  def getAltAlleleWithHighestAlleleCount() : Allele = {
    // optimization: for bi-allelic sites, just return the 1only alt allele
    if ( alleleContext.isBiallelic() ){ alleleContext.getAlternateAllele(0);}
    else
    {
      val alternateAlleles = alleleContext.getAlternateAlleles()

      val alternateAllelesChrCount =  alternateAlleles.map(genotypeContext.getCalledChrCount(_))
      val highestCount = alternateAllelesChrCount.max
      val highestCountLastIndex = alternateAllelesChrCount.lastIndexOf(highestCount)

      alternateAlleles(highestCountLastIndex)
    }
  }


}


object VariantContext
{
  private val WARN_ABOUT_BAD_END = true;
  val NO_LOG10_PERROR = CommonInfo.NO_LOG10_PERROR;
  val PASSES_FILTERS : Set[String] = immutable.HashSet[String]()
  val NO_GENOTYPES = GenotypesContext.NO_GENOTYPES;
  private val NO_VALIDATION = ValidationMode.ValueSet.empty
}

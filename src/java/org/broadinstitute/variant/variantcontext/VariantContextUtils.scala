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

;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;

import org.broadinstitute.variant.utils.GeneralUtils;
import org.broadinstitute.variant.vcf._

import scala.collection._
import org.broad.tribble.TribbleException
import org.apache.commons.jexl2.JexlEngine
import org.apache.commons.jexl2.Expression
import net.sf.samtools.util.Lazy


object VariantContextUtils {

  private val MISSING_KEYS_WARNED_ABOUT = mutable.HashSet[String]()
  private val ASSUME_MISSING_FIELDS_ARE_STRINGS = false;

  /** Use a {@link Lazy} {@link JexlEngine} instance to avoid class-loading issues. (Applications that access this class are otherwise
    * forced to build a {@link JexlEngine} instance, which depends on some apache logging libraries that mightn't be packaged.) */
  final val engine: Lazy[JexlEngine] = new Lazy[JexlEngine](new Lazy.LazyInitializer[JexlEngine] {
    def make: JexlEngine = {
      val jexl: JexlEngine = new JexlEngine
      jexl.setSilent(false)
      jexl.setLenient(false)
      jexl.setDebug(false)
      return jexl
    }
  })

  /**
   * A simple but common wrapper for matching VariantContext objects using JEXL expressions
   */
   class JexlVCMatchExp( val name : String, val exp : Expression)


  /**
   * Method for creating JexlVCMatchExp from input walker arguments names and exps.  These two arrays contain
   * the name associated with each JEXL expression. initializeMatchExps will parse each expression and return
   * a list of JexlVCMatchExp, in order, that correspond to the names and exps.  These are suitable input to
   * match() below.
   *
   * @param names names
   * @param exps  expressions
   * @return list of matches
   */
   def initializeMatchExps( names : Array[String],  exps : Array[String]) : Array[JexlVCMatchExp] = {
    if ( names == null || exps == null ){
      throw new IllegalArgumentException("BUG: neither names nor exps can be null: names " + names + " exps=" + exps )}

    if ( names.length != exps.length ){
      throw new IllegalArgumentException("Inconsistent number of provided filter names and expressions: names=" + names + " exps=" + exps)}


    val map = names.zip(exps).toMap

    VariantContextUtils.initializeMatchExps(map)
  }

  /**
   * Method for creating JexlVCMatchExp from input walker arguments mapping from names to exps.  These two arrays contain
   * the name associated with each JEXL expression. initializeMatchExps will parse each expression and return
   * a list of JexlVCMatchExp, in order, that correspond to the names and exps.  These are suitable input to
   * match() below.
   *
   * @param names_and_exps mapping of names to expressions
   * @return list of matches
   */
  def initializeMatchExps( names_and_exps  :Map[String, String]) : Array[JexlVCMatchExp] = {
    val exps = mutable.ArrayBuffer[JexlVCMatchExp]()

    for ( (name,expStr )  <- names_and_exps ) {

      if ( name == null || expStr == null ) throw new IllegalArgumentException("Cannot create null expressions : " + name +  " " + expStr);
      try {
        val exp = engine.get().createExpression(expStr)
        exps += new JexlVCMatchExp(name, exp)
      } catch{case e: Exception =>  throw new IllegalArgumentException("Argument " + name + "has a bad value. Invalid expression used (" + expStr + "). Please see the JEXL docs for correct syntax.") }
    }

    exps.toArray
  }

  /**
   * Returns true if exp match VC.  See collection<> version for full docs.
   * @param vc    variant context
   * @param exp   expression
   * @return true if there is a match
   */
  def matchExpr( vc : VariantContext,  exp  : JexlVCMatchExp) : Boolean = {
    matchExpr(vc,Array(exp)).get(exp)
  }

  /**
   * Matches each JexlVCMatchExp exp against the data contained in vc, and returns a map from these
   * expressions to true (if they matched) or false (if they didn't).  This the best way to apply JEXL
   * expressions to VariantContext records.  Use initializeMatchExps() to create the list of JexlVCMatchExp
   * expressions.
   *
   * @param vc   variant context
   * @param exps expressions
   * @return true if there is a match
   */
  def matchExpr( vc : VariantContext, exps : Iterable[JexlVCMatchExp]) : JEXLMap = {
     new JEXLMap(scala.collection.JavaConversions.asJavaCollection(exps),vc)

  }


  /**
   * Returns true if exp match VC/g.  See collection<> version for full docs.
   * @param vc   variant context
   * @param g    genotype
   * @param exp   expression
   * @return true if there is a match
   */
  def matchExpr( vc : VariantContext,  g : Genotype, exp : JexlVCMatchExp) : Boolean = {
    matchExpr(vc,g,Array(exp)).get(exp)
  }

  /**
   * Matches each JexlVCMatchExp exp against the data contained in vc/g, and returns a map from these
   * expressions to true (if they matched) or false (if they didn't).  This the best way to apply JEXL
   * expressions to VariantContext records/genotypes.  Use initializeMatchExps() to create the list of JexlVCMatchExp
   * expressions.
   *
   * @param vc   variant context
   * @param g    genotype
   * @param exps expressions
   * @return true if there is a match
   */
  def matchExpr( vc : VariantContext, g : Genotype, exps : Iterable[JexlVCMatchExp]) : JEXLMap =  {
     new JEXLMap(scala.collection.JavaConversions.asJavaCollection(exps),vc,g)
  }



  /**
   * Computes the alternate allele frequency at the provided {@link VariantContext} by dividing its "AN" by its "AC".
   * @param vc The variant whose alternate allele frequency is computed
   * @return The alternate allele frequency in [0, 1]
   * @throws AssertionError When either annotation is missing, or when the compuated frequency is outside the expected range
   */
  def calculateAltAlleleFrequency(vc: VariantContext): Double = {
    if (!vc.hasAttribute(VCFConstants.ALLELE_NUMBER_KEY) || !vc.hasAttribute(VCFConstants.ALLELE_COUNT_KEY)) {
      val AN_KEY = VCFConstants.ALLELE_NUMBER_KEY
      val AC_KEY = VCFConstants.ALLELE_COUNT_KEY
      throw new AssertionError(s"Cannot compute the provided variant's alt allele frequency because it does not have both $AN_KEY and $AC_KEY annotations: $vc")
    }

    val altAlleleCount = vc.getAttributeAsInt(VCFConstants.ALLELE_COUNT_KEY, 0)
    val totalCount = vc.getAttributeAsInt(VCFConstants.ALLELE_NUMBER_KEY, 0)
    val aaf = altAlleleCount / totalCount
    if (aaf > 1 || aaf < 0) {
      throw new AssertionError(s"Expected a minor allele frequency in the range [0, 1], but got $aaf. vc=$vc")
    }
    aaf
  }

  /**
   * Update the attributes of the attributes map given the VariantContext to reflect the
   * proper chromosome-based VCF tags
   *
   * @param vc          the VariantContext
   * @param attributes  the attributes map to populate; must not be null; may contain old values
   * @param removeStaleValues should we remove stale values from the mapping?
   * @return the attributes map provided as input, returned for programming convenience
   */
  def calculateChromosomeCounts(vc: VariantContext, attributes: immutable.Map[String, Any], removeStaleValues: Boolean): Map[String, Any] = {
    calculateChromosomeCounts(vc, attributes, removeStaleValues, mutable.HashSet[String]())
  }

  /**
   * Update the attributes of the attributes map given the VariantContext to reflect the
   * proper chromosome-based VCF tags
   *
   * @param vc          the VariantContext
   * @param attributes  the attributes map to populate; must not be null; may contain old values
   * @param removeStaleValues should we remove stale values from the mapping?
   * @param founderIds - Set of founders Ids to take into account. AF and FC will be calculated over the founders.
   *                   If empty or null, counts are generated for all samples as unrelated individuals
   * @return the attributes map provided as input, returned for programming convenience
   */
  def calculateChromosomeCounts(vc: VariantContext, attributes: immutable.Map[String, Any], removeStaleValues: Boolean, founderIds: Set[String]): immutable.Map[String, Any] = {

    var updatedAttributes = attributes

    val AN = vc.genotypeContext.getCalledChrCount()

    // if everyone is a no-call, remove the old attributes if requested
    if (AN == 0 && removeStaleValues) {
      if (updatedAttributes.contains(VCFConstants.ALLELE_COUNT_KEY))
        updatedAttributes -= VCFConstants.ALLELE_COUNT_KEY
      if (updatedAttributes.contains(VCFConstants.ALLELE_FREQUENCY_KEY))
        updatedAttributes -= VCFConstants.ALLELE_FREQUENCY_KEY
      if (updatedAttributes.contains(VCFConstants.ALLELE_NUMBER_KEY))
        updatedAttributes -= VCFConstants.ALLELE_NUMBER_KEY
      updatedAttributes;
    }
    else {

      if (vc.hasGenotypes) {
        updatedAttributes += (VCFConstants.ALLELE_NUMBER_KEY -> AN)

        // if there are alternate alleles, record the relevant tags
        if (vc.alleleContext.getAlternateAlleles().size > 0) {
          val altAlleles = vc.alleleContext.getAlternateAlleles()

          val alleleCounts = altAlleles.map(a => vc.genotypeContext.getCalledChrCount(a))
          val foundersAlleleCounts = altAlleles.map(a => vc.genotypeContext.getCalledChrCount(a, founderIds))
          val totalFoundersChromosomes = vc.genotypeContext.getCalledChrCount(founderIds).toDouble

          val alleleFreqs = if (AN == 0) {
            foundersAlleleCounts.map(a => 0.0)
          } else {
            foundersAlleleCounts.map(_ / totalFoundersChromosomes)
          }


          updatedAttributes += (VCFConstants.ALLELE_COUNT_KEY -> (if (alleleCounts.size == 1) {
            alleleCounts(0)
          } else {
            alleleCounts
          }))
          updatedAttributes += (VCFConstants.ALLELE_FREQUENCY_KEY -> (if (alleleFreqs.size == 1) {
            alleleFreqs(0)
          } else {
            alleleFreqs
          }))
        } else {
          // if there's no alt AC and AF shouldn't be present
          updatedAttributes -= VCFConstants.ALLELE_COUNT_KEY
          updatedAttributes -= VCFConstants.ALLELE_FREQUENCY_KEY
        }
      }

      updatedAttributes;
    }
  }

  /**
   * Update the attributes of the attributes map in the VariantContextBuilder to reflect the proper
   * chromosome-based VCF tags based on the current VC produced by builder.make()
   *
   * @param builder     the VariantContextBuilder we are updating
   * @param removeStaleValues should we remove stale values from the mapping?
   */
  def calculateChromosomeCounts(builder: VariantContextBuilder, removeStaleValues: Boolean) {
    val vc = builder.make();
    builder.attributes(calculateChromosomeCounts(vc, vc.getAttributes, removeStaleValues, mutable.HashSet[String]()))
  }

  /**
   * Update the attributes of the attributes map in the VariantContextBuilder to reflect the proper
   * chromosome-based VCF tags based on the current VC produced by builder.make()
   *
   * @param builder     the VariantContextBuilder we are updating
   * @param founderIds - Set of founders to take into account. AF and FC will be calculated over the founders only.
   *                   If empty or null, counts are generated for all samples as unrelated individuals
   * @param removeStaleValues should we remove stale values from the mapping?
   */
  def calculateChromosomeCounts(builder: VariantContextBuilder, removeStaleValues: Boolean, founderIds: Set[String]) {
    val vc = builder.make();
    builder.attributes(calculateChromosomeCounts(vc, vc.getAttributes, removeStaleValues, founderIds))
  }

  def getMetaDataForField(header: VCFHeader, field: String): VCFCompoundHeaderLine = {
    var metaData: VCFCompoundHeaderLine = header.getFormatHeaderLine(field);
    if (metaData == null) metaData = header.getInfoHeaderLine(field);
    if (metaData == null) {
      if (ASSUME_MISSING_FIELDS_ARE_STRINGS) {
        if (!MISSING_KEYS_WARNED_ABOUT.contains(field)) {
          MISSING_KEYS_WARNED_ABOUT.add(field)
          if (GeneralUtils.DEBUG_MODE_ENABLED) {
            System.err.println("Field " + field + " missing from VCF header, assuming it is an unbounded string type")
          }
        }
        new VCFInfoHeaderLine(field, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Auto-generated string header for " + field);
      }
      else {
        throw new TribbleException("Fully decoding VariantContext requires header line for all fields, but none was found for " + field);
      }
    }
    metaData;
  }


  /**
   * Returns a newly allocated VC that is the same as VC, but without genotypes
   * @param vc  variant context
   * @return  new VC without genotypes
   */
  @Requires(Array[String]("vc != null"))
  @Ensures(Array[String]("result != null"))
  def sitesOnlyVariantContext(vc: VariantContext): VariantContext = {
    new VariantContextBuilder(vc).noGenotypes().make()
  }

  /**
   * Returns a newly allocated list of VC, where each VC is the same as the input VCs, but without genotypes
   * @param vcs  collection of VCs
   * @return new VCs without genotypes
   */
  @Requires(Array[String]("vcs != null"))
  @Ensures(Array[String]("result != null"))
  def sitesOnlyVariantContexts(vcs: Seq[VariantContext]): Seq[VariantContext] = {
    vcs.map(sitesOnlyVariantContext(_))
  }

  def getSize(vc: VariantContext): Int = {
    vc.getEnd() - vc.getStart() + 1
  }

  def genotypeNames(genotypes: Seq[Genotype]): Set[String] = {
    genotypes.map(_.getSampleName()).toSet[String]
  }

  /**
   * Compute the end position for this VariantContext from the alleles themselves
   *
   * In the trivial case this is a single BP event and end = start (open intervals)
   * In general the end is start + ref length - 1, handling the case where ref length == 0
   * However, if alleles contains a symbolic allele then we use endForSymbolicAllele in all cases
   *
   * @param alleles the list of alleles to consider.  The reference allele must be the first one
   * @param start the known start position of this event
   * @param endForSymbolicAlleles the end position to use if any of the alleles is symbolic.  Can be -1
   *                              if no is expected but will throw an error if one is found
   * @return this builder
   */
  @Requires(Array[String]("! alleles.isEmpty()", "start > 0", "endForSymbolicAlleles == -1 || endForSymbolicAlleles > 0"))
  def computeEndFromAlleles(alleles: Array[Allele], start: Int, endForSymbolicAlleles: Int): Int = {
    val ref = alleles(0)

    if (ref.isNonReference) {
      throw new IllegalStateException("computeEndFromAlleles requires first allele to be reference")
    }

    if (AlleleContext.hasSymbolicAllelesInArray(alleles)) {
      if (endForSymbolicAlleles == -1) {
        throw new IllegalStateException("computeEndFromAlleles found a symbolic allele but no endForSymbolicAlleles was provided")
      }
      endForSymbolicAlleles
    } else {
      start + Math.max(ref.length() - 1, 0)
    }
  }

}

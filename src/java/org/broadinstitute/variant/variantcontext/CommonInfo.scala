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


import org.broadinstitute.variant.vcf.VCFConstants;
import scala.collection._

/**
 * Common utility routines for VariantContext and Genotype
 *
 * @author depristo
 */

/**
 * Created by Wim Spee on 3/2/14.
 */


class CommonInfo(var name : String, _log10PError : Double,  var filters : Set[String],  private val _attributes : immutable.Map[String, Any]) {



  private var attributes : immutable.Map[String, Any] = if(_attributes != null && !_attributes.isEmpty){ _attributes} else{ CommonInfo.NO_ATTRIBUTES}
  private var log10PError = CommonInfo.NO_LOG10_PERROR
  setLog10PError(_log10PError)



  /**
   * @return the name
   */
  def getName =  { name }

  /**
   * Sets the name
   *
   * @param name    the name associated with this information
   */
  def setName( name : String) {
    if ( name == null ){ throw new IllegalArgumentException("Name cannot be null " + this);}
    this.name = name;
  }


  // ---------------------------------------------------------------------------------------------------------
  //
  // Filter
  //
  // ---------------------------------------------------------------------------------------------------------

  def getFiltersMaybeNull : Set[String] = { filters;  }

  def getFilters : Set[String] = {
    if(filters == null){ CommonInfo.NO_FILTERS} else{ filters.toSet }
  }

  def filtersWereApplied = { filters != null }

  def isFiltered : Boolean =  if(filters == null){ false } else { filters.size > 0 }

  def isNotFiltered  : Boolean =  { ! isFiltered }

  def addFilter( filter : String) {
    if ( filters == null ){ filters = mutable.HashSet[String]() }// immutable -> mutable


    if ( filter == null ) throw new IllegalArgumentException("BUG: Attempting to add null filter " + this)
    if ( getFilters.contains(filter) ) throw new IllegalArgumentException("BUG: Attempting to add duplicate filter " + filter + " at " + this)
    filters += filter
  }

  def addFilters( filtersSeq : Seq[String]) {
    if ( filters == null ) throw new IllegalArgumentException("BUG: Attempting to add null filters at" + this);
    filters ++= filtersSeq
  }

  // ---------------------------------------------------------------------------------------------------------
  //
  // Working with log error rates
  //
  // ---------------------------------------------------------------------------------------------------------

  def hasLog10PError = { getLog10PError != CommonInfo.NO_LOG10_PERROR }


  /**
   * @return the -1 * log10-based error estimate
   */
  def getLog10PError = { log10PError }
  def getPhredScaledQual : Double =  { getLog10PError * -10 }

  def setLog10PError( log10PErrorArg : Double) {
    if ( log10PErrorArg > 0 && log10PErrorArg != CommonInfo.NO_LOG10_PERROR)
    {throw new IllegalArgumentException("BUG: log10PError cannot be > 0 : " + log10PErrorArg)}

    if ( log10PErrorArg.isInfinite  )
    {throw new IllegalArgumentException("BUG: log10PError should not be Infinity")}
    if ( log10PErrorArg.isNaN )
    {throw new IllegalArgumentException("BUG: log10PError should not be NaN")}

    this.log10PError = log10PErrorArg
  }

  // ---------------------------------------------------------------------------------------------------------
  //
  // Working with attributes
  //
  // ---------------------------------------------------------------------------------------------------------
  def clearAttributes() { attributes = immutable.HashMap[String, Any]() }

  /**
   * @return the immutable attribute map
   */
  def  getAttributes() : immutable.Map[String, Any] =  { attributes.toMap }

  // todo -- define common attributes as enum

  def setAttributes( map : Map[String, _]) {
    clearAttributes();
    putAttributes(map);
  }

  def putAttribute( key : String,  value : Any) {
    putAttribute(key, value, false);
  }

  def putAttribute( key : String , value : Any,  allowOverwrites : Boolean) {
    if ( ! allowOverwrites && hasAttribute(key) )
    {throw new IllegalStateException("Attempting to overwrite key->value binding: key = " + key + " this = " + this)}

    if ( attributes == CommonInfo.NO_ATTRIBUTES ) // immutable -> mutable
      attributes = immutable.HashMap[String, Any]();

    attributes += (key -> value);
  }

  def removeAttribute( key : String) {
    if ( attributes == CommonInfo.NO_ATTRIBUTES ) // immutable -> mutable
      attributes = immutable.HashMap[String, Any]();
    attributes -= (key);
  }

  def putAttributes( map : Map[String, _]) {
    if ( map != null ) {
      // for efficiency, we can skip the validation if the map is empty
      if ( attributes.size == 0 ) {
        if ( attributes == CommonInfo.NO_ATTRIBUTES ) // immutable -> mutable
          attributes = immutable.HashMap[String, Any]();
        attributes ++= map;
      } else {
        var counter = 0
        val limit = map.size
        val keyIter = map.keys.iterator
        while(counter <  limit)
        {
          val key = keyIter.next
          putAttribute(key , map(key), false)
          counter+=1
        }
      }
    }
  }

  def hasAttribute( key : String) : Boolean = { attributes.contains(key) }

  def getNumAttributes() : Int = { attributes.size }


  /**
   * @param key    the attribute key
   *
   * @return the attribute value for the given key (or null if not set)
   */
  def getAttribute( key : String) : Any = { attributes(key) }

  def getAttribute( key : String, defaultValue : Any) : Any = {
    if ( hasAttribute(key) ){ attributes(key);}
    else{ defaultValue}
  }

  def  getAttributeAsString( key : String,  defaultValue : String) : String =  {
    val x = getAttribute(key);
    if ( x == null ){ defaultValue}
    else
    {
      if ( x.isInstanceOf[String] ){ x.asInstanceOf[String]}
      else{x.toString}
    }
  }



  def getAttributeAsInt( key : String, defaultValue : Int) : Int = {
    val x = getAttribute(key);
    if ( x == null || x == VCFConstants.MISSING_VALUE_v4 ){ defaultValue}
    else
    {
      if ( x.isInstanceOf[Integer] ){x.asInstanceOf[Integer]}
      else{x.asInstanceOf[String].toInt}; // throws an exception if this isn't a string
    }

  }

  def getAttributeAsDouble( key : String, defaultValue : Double) : Double = {
    val x = getAttribute(key);
    if ( x == null ){ defaultValue }
    else
    {
      if ( x.isInstanceOf[Double] ){ x.asInstanceOf[Double];}
      else
      {
        if(x.isInstanceOf[Integer]){ x.asInstanceOf[Integer].toDouble}
        else{x.asInstanceOf[String].toDouble} // throws an exception if this isn't a string}
      }
    }

  }

  def getAttributeAsBoolean( key : String,  defaultValue : Boolean)  : Boolean = {
    val x = getAttribute(key);
    if ( x == null ) { defaultValue}
    else
    {
      if ( x.isInstanceOf[Boolean]){x.asInstanceOf[Boolean]}
      else{x.asInstanceOf[String].toBoolean}; // throws an exception if this isn't a string}


    }

  }

  //    public String getAttributeAsString(String key)      { return (String.valueOf(getExtendedAttribute(key))); } // **NOTE**: will turn a null Object into the String "null"
  //    public int getAttributeAsInt(String key)            { Object x = getExtendedAttribute(key); return x instanceof Integer ? (Integer)x : Integer.valueOf((String)x); }
  //    public double getAttributeAsDouble(String key)      { Object x = getExtendedAttribute(key); return x instanceof Double ? (Double)x : Double.valueOf((String)x); }
  //    public boolean getAttributeAsBoolean(String key)      { Object x = getExtendedAttribute(key); return x instanceof Boolean ? (Boolean)x : Boolean.valueOf((String)x); }
  //    public Integer getAttributeAsIntegerNoException(String key)  { try {return getAttributeAsInt(key);} catch (Exception e) {return null;} }
  //    public Double getAttributeAsDoubleNoException(String key)    { try {return getAttributeAsDouble(key);} catch (Exception e) {return null;} }
  //    public String getAttributeAsStringNoException(String key)    { if (getExtendedAttribute(key) == null) return null; return getAttributeAsString(key); }
  //    public Boolean getAttributeAsBooleanNoException(String key)  { try {return getAttributeAsBoolean(key);} catch (Exception e) {return null;} }
}

object CommonInfo
{
  val NO_LOG10_PERROR : Double = 1.0;
  val NO_FILTERS = immutable.HashSet[String]()
  val NO_ATTRIBUTES = immutable.HashMap[String, Any]()



}
<?xml version="1.0" encoding="UTF-8"?>

<!-- 
     ==================================================================================== 
                               PUBLIC DOMAIN NOTICE
                National Center for Biotechnology Information
 
        This software/database is a "United States Government Work" under the           
        terms of the United States Copyright Act.  It was written as part of                         
        the author's official duties as a United States Government employee and
        thus cannot be copyrighted.  This software/database is freely available
        to the public for use. The National Library of Medicine and the U.S.
        Government have not placed any restriction on its use or reproduction.
        Although all reasonable efforts have been taken to ensure the accuracy
        and reliability of the software and data, the NLM and the U.S.
        Government do not and cannot warrant the performance or results that
        may be obtained by using this software or data. The NLM and the U.S.
        Government disclaim all warranties, express or implied, including
        warranties of performance, merchantability or fitness for any particular
        purpose.

        Please cite the author in any work or product based on this material.

        Author: Jimmy Jin, Anne Kiang
        File Description: xsl for dbGaP data-dictionary
        Date: 6.25.2007
	Updated Date: 12.05.2008
     ==================================================================================== 
-->



<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<xsl:template match="/">
  <html>
  <body>
    <h2>Data dictionary for data table <xsl:value-of select="/data_table/@id"/>.p<xsl:value-of select="/data_table/@participant_set"/>
      <xsl:if test="/data_table/description">: <xsl:value-of select="/data_table/description"/></xsl:if>
</h2>
    <p>phv# = Phenotype Variable ID, v# = Version Number</p>
    <p>
      <xsl:if test="/data_table/unique_key">
	<b>Unique Keys:</b><br/>
	<xsl:for-each select="/data_table/unique_key">
          <xsl:sort select="attribute::phv"/>phv<xsl:value-of select="@phv"/>.v<xsl:value-of select="@version"/>.p<xsl:value-of select="/data_table/@participant_set"/>&#160;&#160;<xsl:value-of select="unique_key"/>
	  <xsl:value-of select="."/><br/>
	</xsl:for-each>
      </xsl:if>
    </p>

    <table border="1" cellpadding="5">
    

      <tr bgcolor="#ffffff">
            <th align="left">Variable ID</th>
	    <th align="left">Variable</th>
	    <th align="left">Description</th>
	    <th align="left">Type</th>
	    <th align="left">Units</th>
	    <th align="left">Logical minimum</th>
	    <th align="left">Logical maximum</th>
	    <th align="left">Coded values</th>
	    <th align="left">Comment</th>
    </tr>
    <xsl:for-each select="data_table/variable">
      
       <xsl:sort select="attribute::name"/>
      <tr>
        <td><a href="http://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/variable.cgi?study_id={/data_table/@study_id}.p{/data_table/@participant_set}&amp;phv={substring(@id,4,8)+0}"><xsl:value-of select="@id"/>.p<xsl:value-of select="/data_table/@participant_set"/></a></td>
	<td><a>
	    <xsl:attribute name="name">	<xsl:value-of select="@id"/> </xsl:attribute>	
	    <xsl:value-of select="name" />
	</a></td>
	<td><xsl:value-of select="description"/></td>
	<td><xsl:value-of select="type"/>&#160;</td>
	<td><xsl:value-of select="unit"/>&#160;</td>
        <td><xsl:value-of select="logical_min"/>&#160;</td>
        <td><xsl:value-of select="logical_max"/>&#160;</td>
        <td>
          <xsl:for-each select="value">
            <xsl:if test="@code != ''"><xsl:value-of select="@code"/>=</xsl:if><xsl:value-of select="."/>
	    <xsl:if test="position() != last()"><br/>
            </xsl:if>
          </xsl:for-each>&#160;
        </td>
	<td><xsl:value-of select="comment"/>&#160;</td>
    </tr>
    </xsl:for-each>
    </table>
    <br/>
    <xsl:if test="/data_table/comments">
      <b>Comments:</b><br/>
      <xsl:for-each select="/data_table/comments">
	<xsl:value-of select="."/><br/>
      </xsl:for-each>
    </xsl:if>
  </body>
  </html>
</xsl:template>

</xsl:stylesheet>

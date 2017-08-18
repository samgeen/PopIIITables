/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/simplex/sx_popiiitables.c
 * \date        7/2017
 * \author      Sam Geen
 * \brief       Read + interpolate over PopIII stellar properties 
 * \details     
 * 
 * 
 * \par Major modifications and contributions:
 * 
 * - DD.MM.YYYY Description
 */

/********************************************
 *         Pop III Properties Tables        *
 ********************************************/

struct PopIIITable
{
    // Number of accretion bins
    int nacc;
    // Number of mass bins
    int nmass;
    // Accretion array
    double* accretion;
    // Mass array
    double* mass;
    // Table of values
    double* values;
};

// Read a Pop III table
void sx_readpopiitable(PopIIITable* table, char* filename)
{
    /*
    Inputs:
    table - structure containing table values to pass out
    filename - name of file to open
    */
    // File pointer
    FILE *file;
    // Dummy integer for reading leftover file record bytes
    int dummy;

    // Open file
    file = fopen(filename,"rb");

    // Read the accretion array with Fortran 4-byte file records
    fread(&table.nacc,4,1,file);
    table.accretion = calloc(sizeof(double),table.nacc)
    fread(table.accretion,sizeof(double),table.nacc,file);
    fread(&dummy,4,1,file);

    // Read the mass array 
    fread(&table.nmass,4,1,file);
    table.mass = calloc(sizeof(double),table.nmass)
    fread(table.mass,sizeof(double),table.nmass,file);
    fread(&dummy,4,1,file);

    // Read the value array
    fread(&dummy,4,1,file);
    table.values = calloc(sizeof(double),table.nacc*table.nmass)
    fread(table.values,sizeof(double),table.nacc*table.nmass,file);
    fread(&dummy,4,1,file);

    // Done
    return;    
}
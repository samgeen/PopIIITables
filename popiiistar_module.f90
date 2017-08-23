! /*!
!  * \copyright   This file is part of the AREPO code developed by Volker Springel.
!  * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
!  * \copyright   and contributing authors.
!  *  
!  * \file        src/simplex/popiiistar_module.f90
!  * \date        8/2017
!  * \author      Sam Geen (STG)
!  * \brief       Read + interpolate over PopIII stellar properties 
!  * \details     
!  * 
!  * 
!  * \par Major modifications and contributions:
!  * 
!  * - 23.08.2017 Created (STG)
!  */

! /********************************************
!  *         Pop III Properties Tables        *
!  ********************************************/

MODULE popiiistar_module

! Gets values from popiii star tables
! Sam Geen, August 2017

! Use lookup table module for accretion, mass values
use lookup_table_module

implicit none

public

! Tables for metal cooling with and without flux
type(lookup_table) :: radius_table
type(lookup_table) :: temperature_table

! Unit conversions to cgs
real(dp),parameter::solar_mass_cgs=1.9891d33
real(dp),parameter::solar_mass_per_year_cgs=6.30321217d25
real(dp),parameter::solar_radius_cgs=6.957d10


CONTAINS
!************************************************************************
! Sets up the tables, and then clears them (e.g. on program exit)
SUBROUTINE setup_popiii_tables(dir)
  use amr_commons
  character(len=128)              ::dir,filename
  filename = 'popiii_grid_radius.dat'
  call setup_table(radius_table, filename)
  filename = 'popiii_grid_T.dat'
  call setup_table(temperature_table, filename)
  if (myid .eq. 1) write(*,*) "Set up Pop III stellar evolution tables"
END SUBROUTINE setup_popiii_tables

SUBROUTINE popiii_lookup(accretion,mass,radius,temperature)
  ! Look up (radius, temperature) of a Pop III star based on its accretion rate and mass
  ! RETURNS
  ! 

  real(dp),intent(in)::accretion
  real(dp),intent(in)::mass
  real(dp),intent(out)::radius
  real(dp),intent(out)::temperature
  real(dp)::asolar,msolar
  
  ! Convert input values from cgs to solar
  asolar = accretion/solar_mass_per_year_cgs
  msolar = mass/solar_mass_cgs
  ! Read radius (output = solar radii)
  call find_value2(radius_table,     asolar,msolar,radius)
  ! Read temperature (output = log10(K))
  call find_value2(temperature_table,asolar,msolar,temperature)
  ! Convert output values to cgs
  radius = radius*solar_radius_cgs
  temperature = 10d0**temperature


END SUBROUTINE ssm_radiation


END MODULE
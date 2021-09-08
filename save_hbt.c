#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "allvars.h"
#include "proto.h"

void create_galaxy_files(int snapnum)
{
  // create output files - snapshot option
  int n, i;
  char buf[1000];

  for (n = 0; n < NOUT; n++)
  {
    sprintf(buf, "%s/%s_%03d", OutputDir, FileNameGalaxies, snapnum);
    if (!(FdGalDumps = fopen(buf, "w+")))
    {
      char sbuf[1000];
      sprintf(buf, "can't open file `%s'\n", buf);
      terminate(sbuf);
    }

    fseek(FdGalDumps[n], 1 * sizeof(int), SEEK_SET); /* skip the space for the header */
  }
}

void close_galaxy_files(void)
{
  int n;

  for (n = 0; n < NOUT; n++)
  {
    fseek(FdGalDumps[n], 0, SEEK_SET);
    myfwrite(&TotGalaxies[n], sizeof(int), 1, FdGalDumps[n]);           //Number of trees
    fclose(FdGalDumps[n]);
  }
}

void save_galaxy_append(int filenr, int snapnum, int i, int n)
{
  struct GALAXY_OUTPUT galaxy_output;

  prepare_galaxy_for_output(n, filenr, snapnum, &HaloGal[i], &galaxy_output);

  myfwrite(&galaxy_output, sizeof(struct GALAXY_OUTPUT), 1, FdGalDumps[n]);

  TotGalaxies[n]++;     //this will be written later
}

void prepare_galaxy_for_output(int n, int filenr, int tree, struct GALAXY *g, struct GALAXY_OUTPUT *o)
{
  int j;

#ifdef OUTPUT_MET_AGE_DISTR
  int kk;
#endif

  /* #ifdef GALAXYTREE */
  /*   long long big; */
  /*   double scalefac, xx, yy, zz; */
  /*   int i, k, lenmax, tmpfirst, next; */
  /* #endif */

#ifndef NO_PROPS_OUTPUTS

  /* #ifdef UPDATETYPETWO */
  /*   if(g->Type == 2) */
  /*     get_coordinates(g->Pos, g->Vel, g->MostBoundID, tree, g->HaloNr, g->SnapNum); */
  /* #endif */

#ifdef NOIRA
  int imet;
#ifdef OUTPUT_MET_DISTR
  int i;
#endif
#endif

  o->Type = g->Type;
#ifndef GALAXYTREE
  long long big = calc_big_db_offset(filenr, tree); //xielizhi
  o->HaloIndex = g->HaloNr;
  o->FirstHaloInFOFgroup = Halo[g->HaloNr].FirstHaloInFOFgroup;

#endif
  o->SnapNum = g->SnapNum;
  o->CentralMvir = get_virial_mass(Halo[g->HaloNr].FirstHaloInFOFgroup);

  for (j = 0; j < 3; j++)
  {
    o->Pos[j] = g->Pos[j];
    o->Vel[j] = g->Vel[j];
    o->GasSpin[j] = g->GasSpin[j];
    o->StellarSpin[j] = g->StellarSpin[j];
  }

  o->GalID = g->GalID;
  o->Len = g->Len;
  o->Mvir = g->Mvir;
  o->Rvir = g->Rvir;
  o->Vvir = g->Vvir;
  o->Vmax = g->Vmax;

  o->ColdGas = g->ColdGas;
  o->H_mol = g->H_mol; //xielizhi
  o->StellarMass = g->StellarMass;
  o->BulgeMass = g->BulgeMass;
  o->HotGas = g->HotGas;
  o->EjectedMass = g->EjectedMass;
  o->BlackHoleMass = g->BlackHoleMass;
  //  o->QuasarModeRate = g->QuasarModeRate;  //xielizhi 20150720 test remember to remove
  //  o->EddingtonRate = g->EddingtonRate;  //xielizhi 20150720 test remember to remove
  //  o->RadioModeRate = g->RadioModeRate; //xielizhi 20150720 test remember to remove

#ifdef NOIRA
  for (imet = 0; imet < NMET; imet++)
  {
    o->NMetalsColdGas[imet] = g->NMetalsColdGas[imet];
    o->NMetalsStellarMass[imet] = g->NMetalsStellarMass[imet];
    o->NMetalsBulgeMass[imet] = g->NMetalsBulgeMass[imet];
    o->NMetalsHotGas[imet] = g->NMetalsHotGas[imet];
    o->NMetalsEjectedMass[imet] = g->NMetalsEjectedMass[imet];
  }
#ifdef NOIRA_STORE_IaRATES
#ifdef NOIRA_STORE_IaRATES_full
  for (j = 0; j < MAXOUTPUTS; j++)
    o->IaRATE = g->SEP.IaRATE[j] * UnitMass_in_g / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS;
#else
  o->IaRATE = g->SEP.IaRATE[n] * UnitMass_in_g / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS;
#endif
#endif
#else
  o->MetalsColdGas = g->MetalsColdGas;
  o->MetalsStellarMass = g->MetalsStellarMass;
  o->MetalsBulgeMass = g->MetalsBulgeMass;
  o->MetalsHotGas = g->MetalsHotGas;
  o->MetalsEjectedMass = g->MetalsEjectedMass;
#endif

    /*  NOTE: in Msun/yr */
#ifdef UseFullSfr
  /* changes in UseFullSfr by Starkenburg, no SfrBulge and adapted to save at every 1,2,4,8 snapshot in accordance with the treefiles */
  o->Mstarsformed = g->Mstarsformed;
  for (j = 0; j < MAXOUTPUTS; j++)
  {
    o->Sfr[j] = g->Sfr[j] * UnitMass_in_g / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS;
    o->SfrBulge[j] = g->SfrBulge[j] * UnitMass_in_g / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS;
  }
#else
  o->Sfr = g->Sfr[n] * UnitMass_in_g / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS;
  o->SfrBulge = g->SfrBulge[n] * UnitMass_in_g / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS;
#endif

#ifdef DO_RPS
  o->RPradius = g->RPradius;
#endif
#ifdef RING //xielizhi
#ifdef OUTPUT_RING
  for (j = 0; j < 21; j++) //rbins
  {
    o->SfrRing[j] = g->SfrRing[j];
    o->Ring[j] = g->Ring[j];
    o->FmolRing[j] = g->FmolRing[j];
#ifdef REALJ
    o->GasRing[j] = g->GasRing[j];
#endif
#ifdef DO_RPS
    o->RPRing[j] = g->RPRing[j];
#endif
  }
#endif
#else
#ifdef OUTPUT_RING
//    o->SfrTemp = g->SfrOut;
#endif
#endif

  o->XrayLum = g->XrayLum;
  o->StellarDiskRadius = g->StellarDiskRadius;
  o->GasDiskRadius = g->GasDiskRadius;
  o->BulgeSize = g->BulgeSize;
  o->CoolingRadius = g->CoolingRadius;
#ifdef MARKTREE
  o->Mhot_rps = g->Mhot_rps;
  o->Mhot_tidal = g->Mhot_tidal;
  o->Mhot_infall = g->Mhot_infall;
  o->Mhot_agn = g->Mhot_agn;
  o->Mhot_merger = g->Mhot_merger;
  o->Mhot_stripcold = g->Mhot_stripcold;
  o->Mhot_reheated = g->Mhot_reheated;
  o->Mhot_ejected = g->Mhot_ejected;
  o->Mhot_reincorporate = g->Mhot_reincorporate;
  o->Rs_sf = g->Rs_sf;
  o->Ms_sf = g->Ms_sf;
  o->Rs_majormerger = g->Rs_majormerger;
  o->Ms_majormerger = g->Ms_majormerger;
  o->Rb_majormerger = g->Rb_majormerger;
  o->Mb_majormerger = g->Mb_majormerger;
  o->Rs_minormerger = g->Rs_minormerger;
  o->Ms_minormerger = g->Ms_minormerger;
  o->Rb_minormerger = g->Rb_minormerger;
  o->Mb_minormerger = g->Mb_minormerger;
  o->Rs_instability = g->Rs_instability; // change of whole disk
  o->Rb_instability = g->Rb_instability; // change of whole disk
  o->Mb_instability = g->Mb_instability; // change of stellar disk
  o->Rg_cooling = g->Rg_cooling;
  o->Mg_cooling = g->Mg_cooling;
  o->Rg_sf = g->Rg_sf;
  o->Rg_recycling = g->Rg_recycling;
  o->Mg_recycling = g->Mg_recycling;
  o->Rg_majormerger = g->Rg_majormerger;
  o->Mg_majormerger = g->Mg_majormerger;
  o->Rg_minormerger = g->Rg_minormerger;
  o->Mg_minormerger = g->Mg_minormerger;
  o->Mg_SNfb = g->Mg_SNfb;
  o->Mg_AGNfb = g->Mg_AGNfb;
  o->Mg_rps = g->Mg_rps;
  o->Rg_rps = g->Rg_rps;
#ifdef NOIRA
  for (imet = 0; imet < NMET; imet++)
  {
    o->MetalsHotGas_stripping[imet] = g->MetalsHotGas_stripping[imet];
    o->MetalsHotGas_reincorporate[imet] = g->MetalsHotGas_reincorporate[imet];
    o->MetalsHotGas_stripcold[imet] = g->MetalsHotGas_stripcold[imet];
    o->MetalsHotGas_reheated[imet] = g->MetalsHotGas_reheated[imet];
    o->MetalsHotGas_ejected[imet] = g->MetalsHotGas_ejected[imet];
    o->MetalsHotGas_AGNfb[imet] = g->MetalsHotGas_AGNfb[imet];
    o->MetalsHotGas_merger[imet] = g->MetalsHotGas_merger[imet];
    o->MetalsStar_sf[imet] = g->MetalsStar_sf[imet];
    o->MetalsStar_majormerger[imet] = g->MetalsStar_majormerger[imet];
    o->MetalsStar_minormerger[imet] = g->MetalsStar_minormerger[imet];
    o->MetalsColdGas_cooling[imet] = g->MetalsColdGas_cooling[imet];
    o->MetalsColdGas_recycling[imet] = g->MetalsColdGas_recycling[imet];
    o->MetalsColdGas_majormerger[imet] = g->MetalsColdGas_majormerger[imet];
    o->MetalsColdGas_minormerger[imet] = g->MetalsColdGas_minormerger[imet];
    o->MetalsColdGas_SNfb[imet] = g->MetalsColdGas_SNfb[imet];
    o->MetalsColdGas_AGNfb[imet] = g->MetalsColdGas_AGNfb[imet];
    o->MetalsColdGas_rps[imet] = g->MetalsColdGas_rps[imet];
  }
#endif
#endif
#endif

#ifdef OUTPUT_REST_MAGS
  /* Luminosities in various bands */
  for (j = 0; j < NMAG; j++)
  {
    o->Mag[j] = lum_to_mag(g->Lum[j][n]);
    o->MagBulge[j] = lum_to_mag(g->LumBulge[j][n]);
    o->MagDust[j] = lum_to_mag(g->LumDust[j][n]);
  }
  if (g->StellarMass > 0.0)
  {
    o->MassWeightAge = g->MassWeightAge[n] / g->StellarMass;
    o->MassWeightAge = o->MassWeightAge / 1000. * UnitTime_in_Megayears / Hubble_h;
  }
  else
  {
    o->MassWeightAge = 0.;
  }
#endif

#ifdef OUTPUT_OBS_MAGS
#ifdef COMPUTE_OBS_MAGS
  /* Luminosities in various bands */
  for (j = 0; j < NMAG; j++)
  {
    o->ObsMag[j] = lum_to_mag(g->ObsLum[j][n]);
    o->ObsMagBulge[j] = lum_to_mag(g->ObsLumBulge[j][n]);
    o->ObsMagDust[j] = lum_to_mag(g->ObsLumDust[j][n]);
    o->dObsMag[j] = lum_to_mag(g->dObsLum[j][n]);
    o->dObsMagBulge[j] = lum_to_mag(g->dObsLumBulge[j][n]);
    o->dObsMagDust[j] = lum_to_mag(g->dObsLumDust[j][n]);
  }
#endif
#endif

#ifdef OUTPUT_MET_AGE_DISTR
  for (j = 0; j < NBINS; j++)
    for (kk = 0; kk < NBINS; kk++)
    {
      o->AgeZstars[j][kk] = g->AgeZstars[j][kk];
      o->AgeZstarsbulge[j][kk] = g->AgeZstarsbulge[j][kk];
    }
#endif

#ifdef NOIRA
#ifdef OUTPUT_MET_DISTR
  for (i = 0; i < NBINS; i++)
  {
    for (j = 0; j < NBINS; j++)
    {
      o->OFe_FeH_stars[i][j] = g->OFe_FeH_stars[i][j];
      o->OFe_FeH_bulgestars[i][j] = g->OFe_FeH_bulgestars[i][j];
    }
    //      o->OH_stars[i] = g-> OH_stars[i];
    //      o->FeH_stars[i] = g-> FeH_stars[i];
  }
#endif
#endif

#ifndef NO_PROPS_OUTPUTS
#ifdef GALAXYTREE
  //o->HaloID = HaloIDs[g->HaloNr].HaloID;
  o->Redshift = ZZ[g->SnapNum];

  int ii = (int)floor(o->Pos[0] * ScaleFactor);
  int jj = (int)floor(o->Pos[1] * ScaleFactor);
  int kk = (int)floor(o->Pos[2] * ScaleFactor);

  o->PeanoKey = peano_hilbert_key(ii, jj, kk, Hashbits);

  o->SubID = calc_big_db_subid_index(g->SnapNum, Halo[g->HaloNr].FileNr, Halo[g->HaloNr].SubhaloIndex);

  int tmpfirst = Halo[g->HaloNr].FirstHaloInFOFgroup;
  int lenmax = 0;
  int next = tmpfirst;
  while (next != -1)
  {
    if (Halo[next].Len > lenmax)
    {
      lenmax = Halo[next].Len;
      tmpfirst = next;
    }
    next = Halo[next].NextHaloInFOFgroup;
  }

  o->MMSubID = calc_big_db_subid_index(g->SnapNum, Halo[tmpfirst].FileNr, Halo[tmpfirst].SubhaloIndex);
#endif
#endif
}

long long calc_big_db_offset(int filenr, int treenr)
{
  long long big;
#ifdef MRII
  // TODO assumes < 10^9 galaxies per tree, is that always correct?
  big = (((filenr * (long long)1000000) + treenr) * (long long)1000000000);
#else
  big = (((filenr * (long long)1000000) + treenr) * (long long)1000000);
#endif
  return big;
}

long long calc_big_db_subid_index(int snapnum, int filenr, int subhaloindex)
{
  long long big;
#ifdef MRII
  big = snapnum * (long long)10000000000000 + filenr * (long long)1000000000 + subhaloindex;
#else
  big = snapnum * (long long)1000000000000 + filenr * (long long)100000000 + subhaloindex;
#endif
  return big;
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>

//#ifdef PARALLEL
#if defined(PARALLEL) || defined(FALSEPARALLEL)
#include <mpi.h>
#endif

#include "allvars.h"
#include "proto.h"

int main(int argc, char **argv)
{
	int filenr, treenr, halonr;
	struct stat filestatus;
	FILE *fd;
	char buf[1000];
	time_t start, current;
	//#ifdef PARALLEL
#if defined(PARALLEL) || defined(FALSEPARALLEL)
	//time_t start, current;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
	MPI_Comm_size(MPI_COMM_WORLD, &NTask);
#else
	NTask = 1;
	ThisTask = 0;
#endif //PARALLEL

	if (argc != 2)
	{
		printf("\n  usage: L-Galaxies <parameterfile>\n\n");
		endrun(0);
	}

	// Reads the parameter file, given as an argument at run time
	read_parameter_file(argv[1]);
	mymalloc_init();

	check_options(); // check compatibility of chosen precompiling options

	// for dust model
	mu_seed = -150;

	init();
	time(&start);
#ifdef PARALLEL
	/* a small delay so that processors dont use the same file */
	do
		time(&current);
	while (difftime(current, start) < 10.0 * ThisTask);
#endif

	double time1, time2;
	time1 = second();
	for (filenr = FirstFile; filenr <= LastFile; filenr++)
	{
		time(&start);

		//if((fd = fopen(buf, "w"))) - vsmart
		//fclose(fd);

#ifdef PARALLEL
		printf("\nTask %d reading file nr %d\n\n", ThisTask, filenr);
#endif
		load_tree_table(filenr);



		FILE *fdg = fopen("treengal.dat", "w");
		int snapnum;
		for (snapnum = FirstSnapNr; snapnum < LastSnapShotNr; snapnum++)
		{

			create_galaxy_files(snapnum);

			loadtree(filenr, snapnum);

			gsl_rng_set(random_generator, filenr * 100000 + snapnum);
			//classify_halo;
			if (TotNGal >= MaxGal)
			{
				AllocValue_MaxGal *= ALLOC_INCREASE_FACTOR;
				MaxGal = AllocValue_MaxGal;
				if (MaxGal < ngal + 1)
					MaxGal = ngal + 1;
				Gal = myrealloc_movable(Gal, sizeof(struct GALAXY) * MaxGal);
			}
			construct_galaxies_hbt(filenr, snapnum);
			/*	for(halonr = 0; halonr<TreeNHalos; halonr++)
			if(Halo[halonr].Rank==0)
				construct_galaxies_hbt(filenr,snapnum,halonr);
	*/
			NHaloGal = 0;
		}

		while (NHaloGal)
			output_galaxy(filenr, snapnum);

			//printf("NumGals/MaxGal = %d/%d  %g\n", NumGals, MaxGals, NumGals/((double)MaxGals));

#ifdef GALAXYTREE
		save_galaxy_tree_finalize(filenr, treenr);
		printf("NGalTree = %d TotGalCount = %d\n", NGalTree, TotGalCount);
		fflush(stdout);
		fprintf(fdg, "%d\n", NGalTree);
		//save_galaxy_tree(filenr, tree); - vsmart
		//#else
		//save_galaxies(filenr, tree);
#endif
		save_galaxy_finalize(filenr, snapnum);

		//free_galaxies_and_tree();
	}

	//finalize_galaxy_file(filenr); - vsmart
#ifdef GALAXYTREE
	close_galaxy_tree_file();
#else
	close_galaxy_files();
#endif

#ifdef NOIRA
	int cc;
	for (cc = 0; cc < NMET; cc++)
		printf("\t\t >>> metals ever produced :: %d %g \n",
			   cc, CHECK_ALL_METALS[cc]);

		// last bug check
		//      fclose(bruttecose);
#else
	printf("\t\t >>> metals ever produced :: %g \n",
		   CHECK_ALL_METALS);
#endif
	printf("\t\t >>> maxcasino :: %g\n", maxcasino);

	time(&current);
	printf("\ndone tree file %d in %lds\n\n", filenr, current - start);
	free_tree_table();
}

time2 = second();
printf("%g seconds\n", time2 - time1);
//  write_Stop();
//#ifdef PARALLEL
#if defined(PARALLEL) || defined(FALSEPARALLEL)
MPI_Finalize();
#endif
return 0;
}

void construct_galaxies_hbt(int filenr, int snapnum)
{
	static int halosdone = 0;
	int prog, fofhalo, Gal_index, rank, centralgal;
	int flag = 0;
	int Nsub, halonr, NGalsnap;
	Gal_index = 0;
	for (halonr = 0; halonr < TreeNHalos[snapnum]; halonr++)
	{
		if (Halo[halonr].Rank == 0)
		{
			HaloAux[halonr].DoneFlag = 1;
			halosdone++;
			Halo[halonr].FirstHaloInFOFgroup = halonr;
			Nsub = find_satellite(halonr, Halo[halonr].HostHaloId);
			fofhalo = Halo[halonr].FirstHaloInFOFgroup;
			int SubsinThisHost[Nsub];
			centralgal = GalId;
			while (fofhalo > -1)
			{
				rank = Halo[fofhalo].Rank;
				SubsinThisHost[rank] = fofhalo;
				HaloAux[halonr].DoneFlag = 1;
				halosdone++;
				Gal_index = join_galaxies_of_progenitors(filenr, snapnum, fofhalo, Gal_index, centralgal);
				fofhalo = Halo[fofhalo].NextHaloInFOFgroup;
			}
			NGalsnap += Nsub;

			HaloAux[fofhalo].DoneFlag = 2;

			if (ngal > 0)
			{
				evolve_galaxies(Halo[halonr].FirstHaloInFOFgroup, centralgal, Nsub, filenr, snapnum)
			}
		}
	}
	TotNGal = NGalsnap;
	if (TotNgal != ngal)
	{
		printf("Something's wrong here, count satellite number");
		exit(190);
	}
}

int join_galaxies_of_progenitors(int filenr, int snapnum, int halonr, int ngal, int centralgal)
{
	int i, j, ngal, prog, first_occupied, centralgal;

	prog = find_progenitor(halonr);
	if (prog == -1 && Halo[halonr].Rank == 0 && Halo[halonr].Len > 0)
	{
		init_galaxy(ngal, halonr);
		ngal++;
	}
	else
	{
		if (Halo[halonr].Len == 1)
			continue;
		Gal[ngal].CentralGal = centralgal;
		Gal[ngal] = HaloGal[prog];
		Gal[ngal].Halonr = halonr;
		Gal[ngal].CoolingRadius = 0.0;
		if (Gal[ngal].Type == 0 || Gal[ngal].Type == 1) /* this deals with the central galaxies of subhalos */
		{
			for (j = 0; j < 3; i++)
			{
				Gal[ngal].Pos[j] = Halo[halonr].Pos[j];
				Gal[ngal].Vel[j] = Halo[halonr].Vel[j];
			}
		}
		else
		{
			for (j = 0; j < 3; i++)
			{
				Gal[ngal].Pos[j] = Halo[halonr].MostBoundPos[j];
				Gal[ngal].Vel[j] = Halo[halonr].MostBoundVel[j];
			}
		}
		Gal[ngal].Len = Halo[halonr].Len;
		Gal[ngal].Vmax = Halo[halonr].VmaxPhysical;

		Gal[ngal].Rvir = Halo[halonr].BoundR200CritComoving;
		Gal[ngal].Mvir = Halo[halonr].BoundM200Crit;

		Gal[ngal].HaloSpin = sqrt(Halo[halonr].Spin[0] * Halo[halonr].Spin[0] +
								  Halo[halonr].Spin[1] * Halo[halonr].Spin[1] +
								  Halo[halonr].Spin[2] * Halo[halonr].Spin[2]);

		if (Halo[halonr].SnapNum == Halo[halonr].SnapshotOfLastIsolation)
		{
			Gal[ngal].InfallVmax = Gal[ngal].Vmax;
			Gal[ngal].InfallMvir = Gal[ngal].Mvir;
#if defined(NOSTRANG) || defined(GRADSTRANG)
			Gal[ngal].InfallHotGas = Gal[ngal].HotGas;
			Gal[ngal].InfallRvir = Gal[ngal].Rvir;
			Gal[ngal].InfallVvir = Gal[ngal].Vvir;
			for (j = 0; j < 3; j++)
				Gal[ngal].InfallSpin[j] = Halo[halonr].Spin[j];

			double m200;
			float a1, a2;
			int snap, jm;
			snap = Gal[ngal].SnapNum;
			m200 = log10(Gal[ngal].Mvir) + 10;
			for (jm = 0; jm < 47; jm++)
			{
				if (m200 < MCZM200[snap][jm])
					break;
			}
			if (jm == 0)
				a1 = 1;
			if (jm == 47)
				a1 = 0;
			if (jm > 0 && jm < 47)
			{
				a1 = (MCZM200[snap][jm] - m200) / (MCZM200[snap][jm] - MCZM200[snap][jm - 1]);
			}
			a2 = 1 - a1;
			Gal[ngal].rhos = a1 * MCZRhos[snap][jm - 1][0] + a2 * MCZRhos[snap][jm][0];
			Gal[ngal].rs = a1 * MCZRhos[snap][jm - 1][1] + a2 * MCZRhos[snap][jm][1];
#endif
		}
		if (DiskRadiusMethod < 2)
		{
			Gal[ngal].GasDiskRadius = get_disk_radius(halonr, ngal);
			Gal[ngal].StellarDiskRadius = Gal[ngal].GasDiskRadius;
		}

		if (halonr == Halo[halonr].FirstHaloInFOFgroup)
			Gal[ngal].Type = 0;
		else if (Halo[halonr].Len > 1)
			Gal[ngal].Type = 1;
		else
			Gal[ngal].Type = 2;

		ngal++;
	}

	return ngal;
}

int find_progenitor(int halonr)
{
	int i, GalID, firstprog;
	firstprog = -1;
	if (Halo[halonr].SnapNum > FirstSnapNr)
	{
		GalID = Halo[halonr].TrackId;

		for (i = 0; i < TotNgal; i++)
		{
			if (HaloGal[i].GalID == GalID)
			{
				firstprog = i;
				break;
			}
		}
	}
	return firstprog;
}

int find_satellite(int firsthalo, int hosthalo);
{
	int j;
	int thishalo = firsthalo;
	int nexthalo;
	int NSub = 0;
	for (j = 0; j < totNhalo; j++)
	{
		if (Halo[j].Rank == 0 || HaloAux[j].Doneflag == 1)
			continue;
		else
		{
			if (Halo[j].HostHaloId == hosthalo)
			{
				Halo[j].FirstHaloInFOFgroup = firsthalo;
				nexthalo = j;
				Halo[thishalo].NextHaloInFOFgroup = nexthalo;
				HaloAux[thishalo].Doneflag = 1;
				thishalo = nexthalo;
				NSub++;
			}
		}
	}
	HaloAux[firsthalo].NGalaxies = NSub;
	return NSub;
}

void evolve_galaxies(int halonr, int centralgal, int ngal, int filenr, int snapnum)
{
	int i, p, nstep, merger_centralgal, n, currenthalo, prevgal, start;
	double infallingGas, coolingGas, deltaT, Zcurr, time, temp;

	deltaT = NumToTime(Gal[centralgal].SnapNum) - NumToTime(Halo[halonr].SnapNum);
	Zcurr = ZZ[snapnum];

	if (Gal[centralgal].Type != 0 || Gal[centralgal].HaloNr != halonr)
	{
		printf("centralgal=%d,type=%d,halonr=%d, %d\n", centralgal, Gal[centralgal].Type, Gal[centralgal].HaloNr, halonr);
		printf("Something wrong here ..... \n");
		exit(54);
	}
	/*  basically compute the diff between the hot gas obtained at the end of the 
     previous snapshot and the one obtained at the beginning of the new snapshot
     using the conservation of the baryons 
     NOTE that here we are assuming that the hot gas can be diluted but it is
     also possible that the metallicity increases because the hot gas 
     decreases. If this is reasonable or not .... who knows?  */

#ifndef BHBH
	if (Gal[centralgal].StellarMass + Gal[centralgal].ColdGas > 0 && Gal[centralgal].BlackHoleMass == 0)
		Gal[centralgal].BlackHoleMass = 50 * pow(Gal[centralgal].Mvir, 1.33) / 3e6;
		//1.33 from Volonteri, Natarajan & Guiltekin 2011
#endif

	infallingGas = infall_recipe(centralgal, ngal, Zcurr, deltaT);

#ifdef NOIRA
	set_metaltab(Halo[halonr].SnapNum / (MAXSNAPS / MAXOUTPUTS) - 1);
#endif

	/* we integrate things forward by using a number of intervals equal to STEPS */
	for (nstep = 0; nstep < STEPS; nstep++)
	{
#if !defined(NOSTRANG) && !defined(GRADSTRANG)
		/*  determine cooling gas given halo properties */
		coolingGas = cooling_recipe(centralgal, deltaT / STEPS);
#endif
#ifdef NOIRA
		if (CHECKMETALS_inGAL(centralgal))
		{
			printf("err:: WARN B\n");
			exit(1001);
		}
#endif

		/* test trees - get rid of type 2 galaxies with zero mass */
		//      for(p = 0; p < ngal; p++)
		//	if(Gal[p].Type == 2 && Gal[p].ColdGas + Gal[p].StellarMass < 1.e-8)
		//	  Gal[p].Type = 3;

		/* Loop over all galaxies in the halo */
		for (p = centralgal; p < centralgal + ngal; p++)
		{

			/* don't treat galaxies that have already merged */
			if (Gal[p].Type == 3)
				continue;

			deltaT = NumToTime(Gal[p].SnapNum) - NumToTime(Halo[halonr].SnapNum);
			time = NumToTime(Gal[p].SnapNum) - (nstep + 0.5) * (deltaT / STEPS);

			if (deltaT < 0)
			{
				printf("the impossibile happened! deltaT < 0!\n");
				exit(1111);
			}
#ifdef NOIRA
			if (CHECKMETALS_inGAL(p))
			{
				printf("err:: WARN B1\n");
				exit(1002);
			}
#endif

			/*  for central galaxy only */
			if (p == centralgal)
			{
				add_infall_to_hot(centralgal, infallingGas / STEPS);
#ifdef CHECKCODE
				check_code(centralgal, "central after infall");
#endif

#if defined(NOSTRANG) || defined(GRADSTRANG)
				Gal[centralgal].InfallHotGas = Gal[centralgal].HotGas;
#endif
#ifdef NOIRA
				if (CHECKMETALS_inGAL(p))
				{
					printf("err:: WARN B2\n");
					exit(1002);
				}
#endif

				if (ReIncorporationFactor > 0.0)
					reincorporate_gas(centralgal, deltaT / STEPS);
#ifdef CHECKCODE
				check_code(centralgal, "central after reincorporate");
#endif

#ifdef NOIRA
				if (CHECKMETALS_inGAL(p))
				{
					printf("err:: WARN C\n");
					exit(1002);
				}
#endif

				if (DiskRadiusMethod == 3)
				{
					Gal[centralgal].GasDiskRadius = get_galics_disk_radius(halonr, centralgal, coolingGas);
					Gal[centralgal].StellarDiskRadius = Gal[centralgal].GasDiskRadius;
				}

#if defined(NOSTRANG) || defined(GRADSTRANG)
				/*  determine cooling gas given halo properties */
				coolingGas = cooling_recipe(centralgal, deltaT / STEPS);
#endif
				cool_gas_onto_galaxy(centralgal, coolingGas);
#ifdef CHECKCODE
				check_code(centralgal, "central after cooling");
#endif

#ifdef NOIRA
				if (CHECKMETALS_inGAL(p))
				{
					printf("err:: WARN D\n");
					exit(1003);
				}
#endif
			}
#if defined(NOSTRANG) || defined(GRADSTRANG)
			else
			{
#ifdef EJECTTOSAT
				if (ReIncorporationFactor > 0.0 && Gal[p].EjectedMass > 0)
					reincorporate_gas(p, deltaT / STEPS);
#ifdef CHECKCODE
				check_code(p, "satellite after reincorporate");
#endif
#endif
				if (Gal[p].HotGas > 0)
				{
					coolingGas = cooling_recipe(p, deltaT / STEPS);
					cool_gas_onto_galaxy(p, coolingGas);
#ifdef CHECKCODE
					check_code(p, "satellite after cooling");
#endif
#ifdef NOIRA
					if (CHECKMETALS_inGAL(p))
					{
						printf("err:: WARN D\n");
						exit(1003);
					}
#endif
				}
			}
#endif
			if (Gal[p].ColdGas > 0) //xielizhi
				starformation_and_feedback(p, centralgal, time, deltaT / STEPS, halonr, nstep, Zcurr);
#ifdef CHECKCODE
			check_code(p, "after star formation");
#endif

#ifdef NOIRA
			if (CHECKMETALS_inGAL(p))
			{
				printf("err:: WARN E\n");
				exit(1004);
			}
#endif
		}

		/* Check for merger events */
		for (p = centralgal; p < centralgal + ngal; p++)
		{
			/* orphan galaxies*/
			if (Gal[p].Len == 1)
			{
				Gal[p].Type = 2;
			}
			if (Gal[p].Type == 2)
			{
				merger_centralgal = Halo[Gal[p].HaloNr].SinkTrackId;
				/* if this is the merger time*/
				if (Halo[Gal[p].HaloNr].SnapshotOfSink == snapnum)
				{
					/* find merger_center*/
					for (q = centralgal; q < centralgal + ngal; p++)
					{
						if (Halo[Gal[q].HaloNr].TrackId == merger_centralgal)
						{
							merger_centralgal = q;
							break;
						}
					}
					deal_with_galaxy_merger(p, merger_centralgal, centralgal, time, deltaT, halonr, nstep, Zcurr);
				}
#ifdef CHECKCODE
				check_code(merger_centralgal, "merger_centralgal after merger");
				check_code(centralgal, "centralgal after merger");
#endif
			}
		}
#ifdef NOIRA
		if (CHECKMETALS_inGAL(-ngal))
		{
			printf("err:: WARN F\n");
			exit(1005);
		}
#endif
		//xielizhi Get_H_mol
		for (p = 0; p < ngal; p++)
		{
			Gal[p].H_mol = Get_H_mol(p);
		}
		//xielizhi Get_H_mol
	}
	/* end move forward in interval STEPS */
#ifndef UseFullSfr
	for (p = 0; p < ngal; p++)
	{
		Gal[p].SfrLast = Gal[p].Sfrnow;
		Gal[p].Sfrnow = 0;
	}
#endif
#ifdef NOIRA
	float deltatime;
	for (p = centralgal; p < centralgal + ngal; p++)
	{
		deltaT = (NumToTime(Gal[p].SnapNum) - NumToTime(Halo[halonr].SnapNum));
		deltatime = deltaT / 1000. * UnitTime_in_Megayears / Hubble_h;

		respread_SEP(deltatime, &Gal[p].SEP, SFHTimes, SFHTimeBins, NSFRBINS);
	}
#endif

	/*  if this is an output snapshot apply the dust model to each galaxy */
	for (n = 0; n < NOUT; n++)
	{
		if (Halo[halonr].SnapNum == ListOutputSnaps[n])
		{
			for (p = centralgal; p < centralgal + ngal; p++)
				//dust_model(p, n);
				tsudy(p, n);
			break;
		}
	}
}

void output_galaxy(int filenr, int snapnum)
{
	int p;
	if（TotNgal >= MaxHaloGal)
	{
		int oldmax = MaxHaloGal;
		AllocValue_MaxHaloGal *= ALLOC_INCREASE_FACTOR;
		MaxHaloGal = AllocValue_MaxHaloGal;
		if (MaxHaloGal < TotNgal + 1)
			MaxHaloGal = TotNgal + 1;
		HaloGal = myrealloc_movable(HaloGal, sizeof(struct GALAXY) * MaxHaloGal);
	}
	for (p = 0; p < TotNgal; p++)
	{
		Gal[p].SnapNum = Halo[halonr].SnapNum;
		HaloGal[p] = Gal[p];
	}
	report_memory_usage(&HighMark, "construct_galaxies");
}

void check_options()
{
#ifdef OUTPUT_OBS_MAGS
#ifndef COMPUTE_OBS_MAGS
	printf("> Error : option OUTPUT_OBS MAGS requires option COMPUTE_OBS_MAGS \n");
	exit(32);
#endif
#endif
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "allvars.h"
#include "proto.h"
#include <hdf5.h>
/*Starkenburg: changed the file to open, giving SnapshotsPerOutput as a variable and calculating the right number of snapshots*/


void load_tree_table(int filenr)
{
  int i, n;
  char buf[1000];
  int snap;
  
  TotHalos = 0;
  for(snap =0; snap<LastSnapShotNr+1; snap++)
  {
    sprintf(buf,"%s/SubSnap_%03d.hdf5",SimulationDir,snap);
    if(!(tree_file = H5Fopen(buf,H5F_ACC_RDONLY,H5P_DEFAULT)))
      {
        char sbuf[1000];
        sprintf(sbuf, "can't open file `%s'\n", buf);
        terminate(sbuf);
      }

      hid_t file_id, dataset_id,halo_tid;
      hid_t dataspace;
      hsize_t size[1];
      herr_t err;
      file_id = H5Fopen(buf, H5F_ACC_RDWR, H5P_DEFAULT);
      dataset_id=H5Dopen(file_id, "/Subhalos",H5P_DEFAULT);
      dataspace = H5Dget_space(dataset_id);
      H5Sget_simple_extent_dims(dataspace, size, NULL);
      err = H5Sclose(dataspace);
      err = H5Dclose(dataset_id);
      if(size[0]>0 && TotHalos ==0)
        FirstSnapNr = snap;
      TotHalos += size[0]; 
      FirstHaloInSnap[snap] = TotHalos - size[0];
      TreeNHalos[snap] = size[0];
  }
  if(FirstHaloInSnap[FirstSnapNr] != 0)
    exit(11);


  /*
  if(Ntrees)
    TreeFirstHalo[0] = 0;
  // Define a variable containing the number you have to jump to get from one firshalo to the next. 
  for(i = 1; i < Ntrees; i++)
    TreeFirstHalo[i] = TreeFirstHalo[i - 1] + TreeNHalos[i - 1];
  */

#ifdef PRELOAD_TREES
  Halo_Data = mymalloc("Halo_Data", sizeof(struct halo_data) * totNHalos);
  myfseek(tree_file, sizeof(int) * (2 + Ntrees), SEEK_SET);
  myfread(Halo_Data, totNHalos, sizeof(struct halo_data), tree_file);
#ifdef PARALLEL
  printf("\nTask %d done loading trees_%d\n", ThisTask, filenr);
#endif
#endif

}


void free_tree_table(void)
{

#ifdef PRELOAD_TREES
  myfree(Halo_Data);
#endif

  myfree(TreeNgals[0]);
  
  myfree(TreeFirstHalo);
  myfree(TreeNHalos);

#ifdef UPDATETYPETWO
  myfree(TreeAuxData);
#endif

  fclose(tree_file);
}


void load_tree(int SnapShotNr, int snapnum)
{
  int i, NHalos_snap;
  hid_t file_id,dataset_id,halo_tid;
  hid_t dataspace;
  hsize_t size[1];
  herr_t err;
  char buf[1000];

  sprintf(buf,"%s/SubSnap_%03d.hdf5",SimulationDir,snapnum);
  if(!(tree_file = H5Fopen(buf,H5F_ACC_RDONLY,H5P_DEFAULT)))
    {
      char sbuf[1000];
      sprintf(sbuf, "can't open file `%s'\n", buf);
      terminate(sbuf);
    }
  file_id = H5Fopen(buf, H5F_ACC_RDONLY, H5P_DEFAULT);
  dataset_id=H5Dopen(file_id, "/Subhalos",H5P_DEFAULT);
  dataspace = H5Dget_space(dataset_id);
  H5Sget_simple_extent_dims(dataspace, size, NULL);
  err = H5Sclose(dataspace);
  err = H5Dclose(dataset_id);
  NHalos_snap = size[0];

    //int *TrackId, *Snapshot;
    //TrackId = (int *)malloc(sizeof(int) * size[0]) ;
    //Snapshot = (int *)malloc(sizeof(int) * size[0]);
    //Halo = malloc(sizeof(struct halo_data) * size[0]);
    //halo_data  Halo[size[0]];

    Halo_Data = mymalloc("Halo_Data", sizeof(struct halo_data) * NHalos_snap);
    halo_tid = H5Tcreate(H5T_COMPOUND, sizeof(halo_data));

    H5Tinsert(halo_tid, "TrackId", HOFFSET(halo_data, TrackId),H5T_NATIVE_INT);
    H5Tinsert(halo_tid, "Nbound", HOFFSET(halo_data, Len),H5T_NATIVE_INT);
    H5Tinsert(halo_tid, "Mbound", HOFFSET(halo_data, Mbound), H5T_NATIVE_FLOAT);
    H5Tinsert(halo_tid, "HostHaloId", HOFFSET(halo_data, HostHaloId),H5T_NATIVE_INT);
    H5Tinsert(halo_tid, "Rank", HOFFSET(halo_data, Rank),H5T_NATIVE_INT);
    H5Tinsert(halo_tid, "LastMaxMass", HOFFSET(halo_data, LastMaxMass),H5T_NATIVE_FLOAT);
    H5Tinsert(halo_tid, "SnapshotOfLastMaxMass", HOFFSET(halo_data, SnapshotOfLastMaxMass),H5T_NATIVE_INT);
    H5Tinsert(halo_tid, "SnapshotOfLastIsolation", HOFFSET(halo_data, SnapshotOfLastIsolation),H5T_NATIVE_INT);
    H5Tinsert(halo_tid, "SnapshotOfSink", HOFFSET(halo_data, SnapshotOfSink), H5T_NATIVE_INT);
    H5Tinsert(halo_tid, "VmaxPhysical", HOFFSET(halo_data, Vmax), H5T_NATIVE_FLOAT);
    H5Tinsert(halo_tid, "LastMaxVmaxPhysical", HOFFSET(halo_data, LastMaxVmax), H5T_NATIVE_FLOAT);
    H5Tinsert(halo_tid, "SnapshotOfLastMaxVmax", HOFFSET(halo_data, SnapshotOfLastMaxVmax), H5T_NATIVE_INT);
    H5Tinsert(halo_tid, "BoundR200CritComoving", HOFFSET(halo_data, R200), H5T_NATIVE_FLOAT);
    H5Tinsert(halo_tid, "BoundM200Crit", HOFFSET(halo_data, M200), H5T_NATIVE_FLOAT);
    H5Tinsert(halo_tid, "SinkTrackId", HOFFSET(halo_data, SinkTrackId),H5T_NATIVE_INT);

    hid_t dyn_arr_type;
    hsize_t     locdims[1];
    locdims[0] = 3;
    dyn_arr_type = H5Tarray_create(H5T_NATIVE_FLOAT, 1, locdims);
    //H5Tinsert(halo_tid, "PhysicalAverageVelocity", HOFFSET(halo_data, PhysicalAverageVelocity)+3*sizeof(float), H5T_NATIVE_FLOAT);
    H5Tinsert(halo_tid, "SpecificAngularMomentum", HOFFSET(halo_data, Spin), dyn_arr_type);
    H5Tinsert(halo_tid, "ComovingAveragePosition", HOFFSET(halo_data, Pos), dyn_arr_type);
    H5Tinsert(halo_tid, "PhysicalAverageVelocity", HOFFSET(halo_data, Vel), dyn_arr_type);
    H5Tinsert(halo_tid, "ComovingMostBoundPosition", HOFFSET(halo_data, MostBoundPos), dyn_arr_type);
    H5Tinsert(halo_tid, "PhysicalMostBoundVelocity", HOFFSET(halo_data, MostBoundVel), dyn_arr_type);

    dataset_id = H5Dopen(file_id,"/Subhalos",H5P_DEFAULT);
    err = H5Dread(dataset_id,halo_tid,H5S_ALL,H5S_ALL,H5P_DEFAULT,Halo);
    err = H5Dclose(dataset_id);
    err = H5Fclose(file_id);
  //Allocate HaloAux and Galaxy structures.
  HaloAux = mymalloc("HaloAux", sizeof(struct halo_aux_data) * NHalos_snap);

  for(i = 0; i < NHalos_snap; i++)
    {
      HaloAux[i].DoneFlag = 0;
      HaloAux[i].HaloFlag = 0;
      HaloAux[i].NGalaxies = 0;
      Halo[i].FirstHaloinFOFgroup = -1;
      Halo[i].NextHaloinFOFgroup = -1;
      Halo[i].FirstProgenitor = -1;
      Halo[i].NextProgenitor = -1;
    }

  if(AllocValue_MaxHaloGal == 0)  // back here later -- Aug-16-2021
    AllocValue_MaxHaloGal = 1 + NHalos_snap / (0.25 * MAXSNAPS);

  if(AllocValue_MaxGal == 0)
    AllocValue_MaxGal = 2000;

  MaxHaloGal = AllocValue_MaxHaloGal;
  NHaloGal = 0;
  HaloGal = mymalloc_movable(&HaloGal, "HaloGal", sizeof(struct GALAXY) * MaxHaloGal);
  HaloGalHeap = mymalloc_movable(&HaloGalHeap, "HaloGalHeap", sizeof(int) * MaxHaloGal);
  for(i = 0; i < MaxHaloGal; i++)
    HaloGalHeap[i] = i;

  MaxGal = AllocValue_MaxGal;
  Gal = mymalloc_movable(&Gal, "Gal", sizeof(struct GALAXY) * MaxGal);

#ifdef GALAXYTREE
  if(AllocValue_MaxGalTree == 0)
    AllocValue_MaxGalTree = 1.5 * TreeNHalos[nr];

  MaxGalTree = AllocValue_MaxGalTree;
  GalTree = mymalloc_movable(&GalTree, "GalTree", sizeof(struct galaxy_tree_data) * MaxGalTree);
#endif

}


void free_galaxies_and_tree(void)
{
#ifdef GALAXYTREE
  myfree(GalTree);
#endif
  myfree(Gal);
  myfree(HaloGalHeap);
  myfree(HaloGal);
  myfree(HaloAux);

#ifndef PRELOAD_TREES
  myfree(Halo);
#endif
}

size_t myfread(void *ptr, size_t size, size_t nmemb, FILE * stream)
{
  size_t nread;

  if(size * nmemb > 0)
    {
      if((nread = fread(ptr, size, nmemb, stream)) != nmemb)
	{
	  if(feof(stream))
	    printf("I/O error (fread) has occured: end of file\n");
	  else
	    printf("I/O error (fread) has occured: %s\n", strerror(errno));
	  fflush(stdout);
	  terminate("read error");
	}
    }
  else
    nread = 0;

  return nread;
}

size_t myfwrite(void *ptr, size_t size, size_t nmemb, FILE * stream)
{
  size_t nwritten;

  if(size * nmemb > 0)
    {
      if((nwritten = fwrite(ptr, size, nmemb, stream)) != nmemb)
	    {
	      printf("I/O error (fwrite) has occured: %s\n", strerror(errno));
	      fflush(stdout);
	      terminate("write error");
	    }
    }
  else
    nwritten = 0;

  return nwritten;
}

int myfseek(FILE * stream, long offset, int whence)
{
  if(fseek(stream, offset, whence))
    {
      printf("I/O error (fseek) has occured: %s\n", strerror(errno));
      fflush(stdout);
      terminate("write error");
    }

  return 0;
}



from parameters import *
# Import shutil to copy files
import shutil
# A class for the LigandMDPrep
# TODO Write a submit all script --> for i in {1..4}; do cd $i && sbatch sub_script.sh && cd .. ;done;

class LigandMDPrep:
    # The constructor for the LigandMDPrep class
    def __init__(self):
        """_summary_
        The constructor for the LigandMDPrep class
        """
        pass
    
    # A function that creates a folder
    def create_folder(self, folder_path):
        """_summary_
        Creates a folder
        
        Args:
            folder_path (_type_): _description_
        
        """
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
            print("Folder created: " + folder_path)
        else:
            print("Folder already exists: " + folder_path)

    def create_system_folders(self, 
                                n_replicas, 
                                folder_path,
                                folder_name,
                                topology_file,
                                coordinates_file):
        """_summary_
        Given a number of replicas creates the following folders architecture:
        .
        ├── folder_name  <- Name of the MD experiment
        │   ├── 1    <- folder for replica 1
        │       ├── in <- input files
        │       └── out <- output files
        │   ├── 2    <- folder for replica 2
        │       ├── in <- input files
        │       └── out <- output files
        │   ├── n    <- folder for replica n
        │       ├── in <- input files
        │       └── out <- output files
        
        Args:
            n_replicas (_type_): _description_
            folder_path (_type_): _description_
            folder_name (_type_): _description_
        
        """
        for i in range(1, n_replicas+1):
            self.create_folder(folder_path + "/" + folder_name + "/" + str(i))
            self.create_folder(folder_path + "/" + folder_name + "/" + str(i) + "/in")
            self.create_folder(folder_path + "/" + folder_name + "/" + str(i) + "/out")
            # Copy topology file to the input folder of the system
            shutil.copy(topology_file, folder_path + "/" + folder_name + "/" + str(i) + "/in/topology.top")
            # Copy coordinates file to the input folder of the system
            shutil.copy(coordinates_file, folder_path + "/" + folder_name + "/" + str(i) + "/in/coordinates.crd")

    def restrained_minim(self,
                            folder_path: str,
                            folder_name: str,
                            replica: int,
                            imin: int=LMDP_RESTRAINED_MINIM_IMIN,
                            maxcyc: int=LMDP_RESTRAINED_MINIM_MAXCYC,
                            ncyc: int=LMDP_RESTRAINED_MINIM_NCYC,
                            ntb: int=LMDP_RESTRAINED_MINIM_NTB,
                            ntr: int=LMDP_RESTRAINED_MINIM_NTR,
                            cut: int=LMDP_RESTRAINED_MINIM_CUT,
                            restraintWt: float=LMDP_RESTRAINED_MINIM_RESTRAINT_WT,
                            restraintMask: str=LMDP_RESTRAINED_MINIM_RESTRAINT_MASK):
        """_summary_
        Creates an input file for a replica
        
        Args:
            folder_path (str): _description_
            folder_name (str): _description_
            replica (int): _description_
            imin (int):     Flag to run minimization. 
                            1, Perform an energy minimization, 
                            0, Run molecular dynamics without any minimization.
            maxcyc (int):   The maximum number of cycles of minimization. 
                            Default is 1 to be consistent with Amber’s other minimization
                            facilities although it may be unrealistically short.
            ncyc (int):     If NTMIN is 1 then the method of minimization will be switched
                            from steepest descent to conjugate gradient after NCYC cycles.
            ntb (int):      This variable controls whether or not periodic boundaries are 
                            imposed on the system during the calculation of non-bonded 
                            interactions. Bonds spanning periodic boundaries are not yet 
                            supported.There is no longer any need to set this variable,
                            since it can be determined from igb and ntp
            ntr (int):      Flag for restraining specified atoms in Cartesian space using
                            a harmonic potential, if ntr = 1. The restrained atoms are 
                            determined by the restraintmask string. The force constant 
                            is restraint_wt.The reference coordinates are read in "restrt"
                            format from the "refc" file.
            cut (int):      This is used to specify the nonbonded cutoff, in Angstroms. 
                            For PME, the cutoff is used to limit direct space sum, 
                            and 8.0 is usually a good value. When igb>0, the cutoff is 
                            used to truncate nonbonded pairs (on an atom-by-atom basis);
                            here a larger value than the default is generally required.
            restraintWt (float): The weight (kcal · mol−1 · Å−2) for Cartesian restraints 
                            when ntr = 1. The restraint is of the form k(∆x)2, where k is 
                            the force constant of value given by this variable, and ∆x is 
                            the difference between one of the Cartesian coordinates of a 
                            restrained atom and its reference position. There is a term 
                            like this for each Cartesian coordinate of each restrainted 
                            atom. Note that this variable does not have anything to do 
                            with NMR restraints, and there is no way to have restraint_wt
                            depend upon the time step.
            restraintMask (str): String that specifies the restrained atoms when ntr = 1.
        """
        with open(folder_path + "/" + folder_name + "/" + str(replica) + "/in/restrained_minim.in", "w") as f:
            f.write(f"""
                     &cntrl
                      imin   = {imin},
                      maxcyc = {maxcyc},
                      ncyc   = {ncyc},
                      ntb    = {ntb},
                      ntr    = {ntr},
                      cut    = {cut},
                      restraint_wt={restraintWt},
                      restraintmask='{restraintMask}',
                     /
                    END
                    """.replace("                    ", ""))
            f.close()

    def minim(self,
                folder_path: str,
                folder_name: str,
                replica: int,
                imin: int=LMDP_MINIM_IMIN,
                maxcyc: int=LMDP_MINIM_MAXCYC,
                ncyc: int=LMDP_MINIM_NCYC,
                ntb: int=LMDP_MINIM_NTB,
                ntr: int=LMDP_MINIM_NTR,
                cut: int=LMDP_MINIM_CUT):
        """_summary_
        Creates an input file for a replica
        
        Args:
            folder_path (str): _description_
            folder_name (str): _description_
            replica (int): _description_
            imin (int):     Flag to run minimization. 
                            1, Perform an energy minimization, 
                            0, Run molecular dynamics without any minimization.
            maxcyc (int):   The maximum number of cycles of minimization. 
                            Default is 1 to be consistent with Amber’s other minimization
                            facilities although it may be unrealistically short.
            ncyc (int):     If NTMIN is 1 then the method of minimization will be switched
                            from steepest descent to conjugate gradient after NCYC cycles.
            ntb (int):      This variable controls whether or not periodic boundaries are 
                            imposed on the system during the calculation of non-bonded 
                            interactions. Bonds spanning periodic boundaries are not yet 
                            supported.There is no longer any need to set this variable,
                            since it can be determined from igb and ntp
            ntr (int):      Flag for restraining specified atoms in Cartesian space using
                            a harmonic potential, if ntr = 1. The restrained atoms are 
                            determined by the restraintmask string. The force constant 
                            is restraint_wt.The reference coordinates are read in "restrt"
                            format from the "refc" file.
            cut (int):      This is used to specify the nonbonded cutoff, in Angstroms. 
                            For PME, the cutoff is used to limit direct space sum, 
                            and 8.0 is usually a good value. When igb>0, the cutoff is 
                            used to truncate nonbonded pairs (on an atom-by-atom basis);
                            here a larger value than the default is generally required.
        """
        with open(folder_path + "/" + folder_name + "/" + str(replica) + "/in/minim.in", "w") as f:
            f.write(f"""
                     &cntrl
                      imin   = {imin},
                      maxcyc = {maxcyc},
                      ncyc   = {ncyc},
                      ntb    = {ntb},
                      ntr    = {ntr},
                      cut    = {cut},
                     /
                    END
                     """.replace("                    ", ""))

    def heating1(self,
                    folder_path: str, 
                    folder_name: str,
                    replica: int,
                    imin: int=LMDP_HEATING1_IMIN,
                    irest: int=LMDP_HEATING1_IREST,
                    ntx: int=LMDP_HEATING1_NTX,
                    ntb: int=LMDP_HEATING1_NTB,
                    cut: int=LMDP_HEATING1_CUT,
                    ntr: int=LMDP_HEATING1_NTR,
                    ntc: int=LMDP_HEATING1_NTC,
                    ntf: int=LMDP_HEATING1_NTF,
                    tempi: float=LMDP_HEATING1_TEMPI,
                    temp0: float=LMDP_HEATING1_TEMP0,
                    ntt: int=LMDP_HEATING1_NTT,
                    gamma_ln: float=LMDP_HEATING1_GAMMA_LN,
                    nstlim: int=LMDP_HEATING1_NSTLIM,
                    dt: float=LMDP_HEATING1_DT,
                    ntpr: int=LMDP_HEATING1_NTPR,
                    ntwx: int=LMDP_HEATING1_NTWX,
                    ntwr: int=LMDP_HEATING1_NTWR,
                    restraintmask: str=LMDP_HEATING1_RESTRAINTMASK,
                    restraintWt: int=LMDP_HEATING1_RESTRAINT_WT,
                    ig: int=LMDP_HEATING1_IG):
        """_summary_
        Creates an input file for a replica
        
        Args:
            folder_path (str): _description_
            folder_name (str): _description_
            replica (int): _description_
            imin (int):     Flag to run minimization. 
                            1, Perform an energy minimization, 
                            0, Run molecular dynamics without any minimization.
            irest (int):    Flag to restart from a previous run. 
                            0, Do not restart from a previous run. 
                            1, Restart from a previous run.
            ntx (int):      __description__
            ntb (int):      __description__
            cut (int):      __description__
            ntr (int):      __description__
            ntc (int):      __description__
            ntf (int):      __description__
            tempi (float):  __description__
            temp0 (float):  __description__
            ntt (int):      __description__
            gamma_ln (float): __description__
            nstlim (int):   __description__
            dt (float):     __description__
            ntpr (int):     __description__
            ntwx (int):     __description__
            ntwr (int):     __description__
            restraintmask (str): __description__
            restraintWt (int): __description__
            ig (int):       __description__
        """
        with open(folder_path + "/" + folder_name + "/" + str(replica) + "/in/heating1.in", "w") as f:
            f.write(f"""&cntrl
                      imin   = {imin},
                      irest  = {irest},
                      ntx    = {ntx},
                      ntb    = {ntb},
                      cut    = {cut},
                      ntr    = {ntr},
                      ntc    = {ntc},
                      ntf    = {ntf},
                      tempi  = {tempi},
                      temp0  = {temp0},
                      ntt    = {ntt},
                      gamma_ln = {gamma_ln},
                      nstlim = {nstlim},
                      dt     = {dt},
                      ntpr   = {ntpr},
                      ntwx   = {ntwx},
                      ntwr   = {ntwr},
                      restraintmask = '{restraintmask}',
                      restraint_wt = {restraintWt},
                      ig     = {ig},
                     /
                     """.replace("                    ", ""))
    
    def heating2_0(self, 
                    folder_path: str,
                    folder_name: str,
                    replica: int,
                    imin: int=LMDP_HEATING2_0_IMIN,
                    irest: int=LMDP_HEATING2_0_IREST,
                    ntx: int=LMDP_HEATING2_0_NTX,
                    ntb: int=LMDP_HEATING2_0_NTB,
                    pres0: float=LMDP_HEATING2_0_PRES0,
                    ntp: int=LMDP_HEATING2_0_NTP,
                    taup: float=LMDP_HEATING2_0_TAUP,
                    cut: float=LMDP_HEATING2_0_CUT,
                    ntr: int=LMDP_HEATING2_0_NTR,
                    ntc: int=LMDP_HEATING2_0_NTC,
                    ntf: int=LMDP_HEATING2_0_NTF,
                    tempi: float=LMDP_HEATING2_0_TEMPI,
                    temp0: float=LMDP_HEATING2_0_TEMP0,
                    ntt: int=LMDP_HEATING2_0_NTT,
                    gamma_ln: float=LMDP_HEATING2_0_GAMMA_LN,
                    nstlim: int=LMDP_HEATING2_0_NSTLIM,
                    dt: float=LMDP_HEATING2_0_DT,
                    ntpr: int=LMDP_HEATING2_0_NTPR,
                    ntwx: int=LMDP_HEATING2_0_NTWX,
                    ntwr: int=LMDP_HEATING2_0_NTWR,
                    restraintmask: str=LMDP_HEATING2_0_RESTRAINTMASK,
                    restraintWt: int=LMDP_HEATING2_0_RESTRAINT_WT,
                    ig: int=LMDP_HEATING2_0_IG):
        """_summary_
        Creates an input file for a replica
        
        Args:
            folder_path (str): _description_
            folder_name (str): _description_
            replica (int): _description_
            imin (int):     Flag to run minimization. 
                            1, Perform an energy minimization, 
                            0, Run molecular dynamics without any minimization.
            irest (int):    Flag to restart from a previous run. 
                            0, Do not restart from a previous run. 
                            1, Restart from a previous run.
            ntx (int):      __description__
            ntb (int):      __description__
            pres0 (float):  __description__
            ntp (int):      __description__
            taup (float):   __description__
            cut (int):      __description__
            ntr (int):      __description__
            ntc (int):      __description__
            ntf (int):      __description__
            tempi (float):  __description__
            temp0 (float):  __description__
            ntt (int):      __description__
            gamma_ln (float): __description__
            nstlim (int):   __description__
            dt (float):     __description__
            ntpr (int):     __description__
            ntwx (int):     __description__
            ntwr (int):     __description__
            restraintmask (str): __description__
            restraintWt (int): __description__
            ig (int):       __description__
        """
        with open(folder_path + "/" + folder_name + "/" + str(replica) + "/in/heating2.0.in", "w") as f:
            f.write(f"""&cntrl
                      imin   = {imin},
                      irest  = {irest},
                      ntx    = {ntx},
                      ntb    = {ntb},
                      pres0  = {pres0},
                      ntp    = {ntp},
                      taup   = {taup},
                      cut    = {cut},
                      ntr    = {ntr},
                      ntc    = {ntc},
                      ntf    = {ntf},
                      tempi  = {tempi},
                      temp0  = {temp0},
                      ntt    = {ntt},
                      gamma_ln = {gamma_ln},
                      nstlim = {nstlim},
                      dt     = {dt},
                      ntpr   = {ntpr},
                      ntwx   = {ntwx},
                      ntwr   = {ntwr},
                      restraintmask = '{restraintmask}',
                      restraint_wt = {restraintWt},
                      ig     = {ig},
                     /
                     """.replace("                    ", ""))

    def heating2_1(self, 
                    folder_path: str,
                    folder_name: str,
                    replica: int,
                    imin: int=LMDP_HEATING2_1_IMIN,
                    irest: int=LMDP_HEATING2_1_IREST,
                    ntx: int=LMDP_HEATING2_1_NTX,
                    ntb: int=LMDP_HEATING2_1_NTB,
                    pres0: float=LMDP_HEATING2_1_PRES0,
                    ntp: int=LMDP_HEATING2_1_NTP,
                    taup: float=LMDP_HEATING2_1_TAUP,
                    cut: int=LMDP_HEATING2_1_CUT,
                    ntr: int=LMDP_HEATING2_1_NTR,
                    ntc: int=LMDP_HEATING2_1_NTC,
                    ntf: int=LMDP_HEATING2_1_NTF,
                    tempi: float=LMDP_HEATING2_1_TEMPI,
                    temp0: float=LMDP_HEATING2_1_TEMP0,
                    ntt: int=LMDP_HEATING2_1_NTT,
                    gamma_ln: float=LMDP_HEATING2_1_GAMMA_LN,
                    nstlim: int=LMDP_HEATING2_1_NSTLIM,
                    dt: float=LMDP_HEATING2_1_DT,
                    ntpr: int=LMDP_HEATING2_1_NTPR,
                    ntwx: int=LMDP_HEATING2_1_NTWX,
                    ntwr: int=LMDP_HEATING2_1_NTWR,
                    restraintmask: str=LMDP_HEATING2_1_RESTRAINTMASK,
                    restraintWt: int=LMDP_HEATING2_1_RESTRAINT_WT,
                    ig: int=LMDP_HEATING2_1_IG):
        """
        Save heating2.1.in template input file in the input folder of the system, using f-strings, for each replica.
        Args:
            folder_path (str): __description__
            folder_name (str): __description__
            replica (int): __description__
            imin (int): __description__
            irest (int): __description__
            ntx (int): __description__
            ntb (int): __description__
            pres0 (float): __description__
            ntp (int): __description__
            taup (float): __description__
            cut (int): __description__
            ntr (int): __description__
            ntc (int): __description__
            ntf (int): __description__
            tempi (float): __description__
            temp0 (float): __description__
            ntt (int): __description__
            gamma_ln (float): __description__
            nstlim (int): __description__
            dt (float): __description__
            ntpr (int): __description__
            ntwx (int): __description__
            ntwr (int): __description__
            restraintmask (str): __description__
            restraintWt (int): __description__
            ig (int): __description__
        """
        with open(folder_path + "/" + folder_name + "/" + str(replica) + "/in/heating2.1.in", "w") as f:
            f.write(f"""&cntrl
                      imin   = {imin},
                      irest  = {irest},
                      ntx    = {ntx},
                      ntb    = {ntb},
                      pres0  = {pres0},
                      ntp    = {ntp},
                      taup   = {taup},
                      cut    = {cut},
                      ntr    = {ntr},
                      ntc    = {ntc},
                      ntf    = {ntf},
                      tempi  = {tempi},
                      temp0  = {temp0},
                      ntt    = {ntt},
                      gamma_ln = {gamma_ln},
                      nstlim = {nstlim},
                      dt     = {dt},
                      ntpr   = {ntpr},
                      ntwx   = {ntwx},
                      ntwr   = {ntwr},
                      restraintmask = '{restraintmask}',
                      restraint_wt = {restraintWt},
                      ig     = {ig},
                     /
                     """.replace("                    ", ""))

    def heating3(self,
                    folder_path: str,
                    folder_name: str,
                    replica: int,
                    imin: int=LMDP_HEATING3_IMIN,
                    irest: int=LMDP_HEATING3_IREST,
                    ntx: int=LMDP_HEATING3_NTX,
                    ntb: int=LMDP_HEATING3_NTB,
                    pres0: float=LMDP_HEATING3_PRES0,
                    ntp: int=LMDP_HEATING3_NTP,
                    taup: float=LMDP_HEATING3_TAUP,
                    cut: int=LMDP_HEATING3_CUT,
                    ntr: int=LMDP_HEATING3_NTR,
                    ntc: int=LMDP_HEATING3_NTC,
                    ntf: int=LMDP_HEATING3_NTF,
                    tempi: float=LMDP_HEATING3_TEMPI,
                    temp0: float=LMDP_HEATING3_TEMP0,
                    ntt: int=LMDP_HEATING3_NTT,
                    gamma_ln: float=LMDP_HEATING3_GAMMA_LN,
                    nstlim: int=LMDP_HEATING3_NSTLIM,
                    dt: float=LMDP_HEATING3_DT,
                    ntpr: int=LMDP_HEATING3_NTPR,
                    ntwx: int=LMDP_HEATING3_NTWX,
                    ntwr: int=LMDP_HEATING3_NTWR,
                    ig: int=LMDP_HEATING3_IG):
        """
        Save heating3.in template input file in the input folder of the system, use f-strings, for each replica.
        Args:
            folder_path (str): __description__
            folder_name (str): __description__
            replica (int): __description__
            imin (int): __description__
            irest (int): __description__
            ntx (int): __description__
            ntb (int): __description__
            pres0 (float): __description__
            ntp (int): __description__
            taup (float): __description__
            cut (int): __description__
            ntr (int): __description__
            ntc (int): __description__
            ntf (int): __description__
            tempi (float): __description__
            temp0 (float): __description__
            ntt (int): __description__
            gamma_ln (float): __description__
            nstlim (int): __description__
            dt (float): __description__
            ntpr (int): __description__
            ntwx (int): __description__
            ntwr (int): __description__
            ig (int): __description__
        """
        with open(folder_path + "/" + folder_name + "/" + str(replica) + "/in/heating3.in", "w") as f:
            f.write(f"""&cntrl
                      imin   = {imin},
                      irest  = {irest},
                      ntx    = {ntx},
                      ntb    = {ntb},
                      pres0  = {pres0},
                      ntp    = {ntp},
                      taup   = {taup},
                      cut    = {cut},
                      ntr    = {ntr},
                      ntc    = {ntc},
                      ntf    = {ntf},
                      tempi  = {tempi},
                      temp0  = {temp0},
                      ntt    = {ntt},
                      gamma_ln = {gamma_ln},
                      nstlim = {nstlim},
                      dt     = {dt},
                      ntpr   = {ntpr},
                      ntwx   = {ntwx},
                      ntwr   = {ntwr},
                      ig     = {ig},
                     /
                      """.replace("                    ", ""))

    def md(self,
            folder_path: str,
            folder_name: str,
            replica: int,
            imin: int=LMDP_MD_IMIN,
            irest: int=LMDP_MD_IREST,
            ntx: int=LMDP_MD_NTX,
            cut: float=LMDP_MD_CUT,
            ntr: int=LMDP_MD_NTR,
            ntb: int=LMDP_MD_NTB,
            ntp: int=LMDP_MD_NTP,
            taup: float=LMDP_MD_TAUP,
            ntc: int=LMDP_MD_NTC,
            ntf: int=LMDP_MD_NTF,
            temp0: float=LMDP_MD_TEMP0,
            ntt: int=LMDP_MD_NTT,
            gamma_ln: float=LMDP_MD_GAMMA_LN,
            nstlim: int=LMDP_MD_NSTLIM,
            dt: float=LMDP_MD_DT,
            ntpr: int=LMDP_MD_NTPR,
            ntwx: int=LMDP_MD_NTWX,
            ntwr: int=LMDP_MD_NTWR,
            ig: int=LMDP_MD_IG):
        """
        Save md.in template input file in the input folder of the system, use f-strings, for each replica.
        Args:
            folder_path (str): __description__
            folder_name (str): __description__
            replica (int): __description__
            imin (int): __description__
            irest (int): __description__
            ntx (int): __description__
            cut (float): __description__
            ntr (int): __description__
            ntb (int): __description__
            ntp (int): __description__
            taup (float): __description__
            ntc (int): __description__
            ntf (int): __description__
            temp0 (float): __description__
            ntt (int ): __description__
            gamma_ln (float): __description__
            nstlim (int): __description__
            dt (float): __description__
            ntpr (int): __description__
            ntwx (int): __description__
            ntwr (int): __description__
            ig (int): __description__
        """
        with open(folder_path + "/" + folder_name + "/" + str(replica) + "/in/md.in", "w") as f:
            f.write(f"""&cntrl
                      imin   = {imin},
                      irest  = {irest},
                      ntx    = {ntx},
                      cut    = {cut},
                      ntr    = {ntr},
                      ntb    = {ntb},
                      ntp    = {ntp},
                      taup   = {taup},
                      ntc    = {ntc},
                      ntf    = {ntf},
                      temp0  = {temp0},
                      ntt    = {ntt},
                      gamma_ln = {gamma_ln},
                      nstlim = {nstlim},
                      dt     = {dt},
                      ntpr   = {ntpr},
                      ntwx   = {ntwx},
                      ntwr   = {ntwr},
                      ig     = {ig},
                     /
                      """.replace("                    ", ""))

    def slurm(self, 
                            folder_path: str,
                            folder_name: str,
                            replica: int,
                            partition: str,
                            time: str,
                            gpu: int,
                            email: str,
                            mod: str='Amber/20-GCC-8.3.0-CUDA-10.1.243'):
        """
        Create slurm file for each replica.
        Args:
            folder_path (str): __description__
            folder_name (str): __description__
            replica (int): __description__
            partition (str): __description__
            time (str): __description__
            gpu (int): __description__
            email (str): __description__
        """
        with open(folder_path + "/" + folder_name + "/" + str(replica) + "/sub_script.sh" , "w") as f:
            f.write(f"""#!/bin/bash
                    #SBATCH --job-name={replica}{folder_name}
                    #SBATCH --partition={partition}
                    #SBATCH --time={time}
                    #SBATCH --gres=gpu:{gpu}:1
                    #SBATCH --nodes=1
                    #SBATCH --ntasks-per-node=1
                    #SBATCH --mail-user={email}
                    #SBATCH --mail-type=ALL

                    module load {mod}

                    echo "rest_minim" > log
                    pmemd.cuda -O -i ./in/restrained_minim.in \\
                                  -o ./out/{folder_name}_restrained.out \\
                                  -p ./in/topology.top -c ./in/coordinates.crd \\
                                  -r ./out/{folder_name}_restrained.rst \\
                                  -ref ./in/coordinates.crd &&
                    echo "minim" >> log
                    pmemd.cuda -O \\
                               -i ./in/minim.in \\
                               -o ./out/{folder_name}_restrained2.out \\
                               -p ./in/topology.top \\
                               -c ./out/{folder_name}_restrained.rst \\
                               -r ./out/{folder_name}_restrained2.rst \\
                               -ref ./in/coordinates.crd &&
                    echo "heat1" >> log
                    pmemd.cuda -O \\
                               -i ./in/heating1.in \\
                               -o ./out/{folder_name}_heat1.out \\
                               -p ./in/topology.top \\
                               -c ./out/{folder_name}_restrained2.rst \\
                               -r ./out/{folder_name}_heat1.rst \\
                               -ref ./out/{folder_name}_restrained2.rst \\
                               -x ./out/heat1.nc &&
                    echo "heat2.0" >> log
                    pmemd.cuda -O \\
                               -i ./in/heating2.0.in \\
                               -o ./out/{folder_name}_heat2.0.out \\
                               -p ./in/topology.top \\
                               -c ./out/{folder_name}_heat1.rst \\
                               -r ./out/{folder_name}_heat2.0.rst \\
                               -x ./out/heat2.0.nc &&
                    echo "heat2.1" >> log
                    pmemd.cuda -O \\
                               -i ./in/heating2.1.in \\
                               -o ./out/{folder_name}_heat2.1.out \\
                               -p ./in/topology.top \\
                               -c ./out/{folder_name}_heat2.0.rst \\
                               -r ./out/{folder_name}_heat2.1.rst \\
                               -x ./out/heat2.1.nc &&
                    echo "heat3" >> log
                    pmemd.cuda -O \\
                               -i ./in/heating3.in \\
                               -o ./out/{folder_name}_heat3.out \\
                               -p ./in/topology.top \\
                               -c ./out/{folder_name}_heat2.1.rst \\
                               -r ./out/{folder_name}_heat3.rst \\
                               -x ./out/heat3.nc &&
                    echo "md" >> log
                    pmemd.cuda -O \\
                               -i ./in/md.in \\
                               -o ./out/{folder_name}_md.out \\
                               -p ./in/topology.top \\
                               -c ./out/{folder_name}_heat3.rst \\
                               -r ./out/{folder_name}_md.rst \\
                               -x ./out/md.nc""".replace("                    ", ""))

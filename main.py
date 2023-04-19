import os,math,random,shutil,time,queue
from multiprocessing import Process,Queue
import thread
import results as postprocess
# source /home/artem/cp2k/tools/toolchain/install/setup


TOP = 'atropoisomer.prmtop'
PDB = 'atropoisomer.pdb'
PROJECT_NAME = 'atropo'
XYZ = 'atropoisomer-1.xyz'


# os.chdir('/home/artem/HREMD_algorithm/cp2k/dinuclear-Cu2L2/')
os.chdir('/home/artem/HREMD_algorithm/cp2k/')

# TOP = 'Cu2L2.prmtop'
# PDB = 'Cu2L2.pdb'
# PROJECT_NAME = 'Cu2L2_coord'
# XYZ = 'Cu2L2.xyz'
# RMSD_REF = '/home/artem/HREMD_algorithm/cp2k/dinuclear-Cu2L2/Cu2L2-ref.xyz'
RMSD_REF = '/home/artem/HREMD_algorithm/cp2k/rmsd/atropo2-cp2kopt.xyz'

ATOMS_EXT = "1..66"
#ATOMS_EXT = "1..118"
# ATOMS_EXT = "1..13 15..17 19 21 23 28 30 32 37..39 41 43 45 48..60 62..64 66 68 70 75 77 79 84..86 88 90 92 95..96"

alpha_vs = 0.01 # default: 0.15
delta_vs = 0.10 # default: 0.10
steps_vs = 1000 # default: 40

md_steps = 2000

TEMPERATURE = 700 # default: 298


PARM = {
    "GLOBAL": [
        {
            "PRINT_LEVEL": [
                "MEDIUM"
            ],
            "PROJECT": [
            f"{PROJECT_NAME}"
            ],
            "RUN_TYPE": [
                "MD"
            ],
            "SEED": [
                "1642"
            ]
        }
    ],
    "FORCE_EVAL": [
        {
            "EXTERNAL_POTENTIAL": [
                {
                    "ATOMS_LIST": [
                        f"{ATOMS_EXT}"
                    ],
                    "FUNCTION": [
                        "0.00000000001*(Z^2)^4"
                    ]
                },
                {
                    "ATOMS_LIST": [
                        f"{ATOMS_EXT}"
                    ],
                    "FUNCTION": [
                        "0.00000000001*(X^2)^4"
                    ]
                },
                {
                    "ATOMS_LIST": [
                        f"{ATOMS_EXT}"
                    ],
                    "FUNCTION": [
                        "0.00000000001*(Y^2)^4"
                    ]
                }
            ],
            "METHOD": [
                "FIST"
            ],
            "MM": [
                {
                    "FORCEFIELD": [
                        {
                            "DO_NONBONDED": [
                                "T"
                            ],
                            "ALL_EI_SCALE": [
                                "1.0"
                            ],
                            "ALL_VDW_SCALE": [
                                "1.0"
                            ],
                            "EI_SCALE14": [
                                "1.0"
                            ],
                            "VDW_SCALE14": [
                                "1.0"
                            ],
                            "FORCE_SCALE": [
                                "1.0"
                            ],
                            "PARMTYPE": [
                                "AMBER"
                            ],
                            "PARM_FILE_NAME": [
                                f"{TOP}"
                            ],
                            "SHIFT_CUTOFF": [
                                "F"
                            ],
                            "SPLINE": [
                                {
                                    "EMAX_SPLINE": [
                                        "1.0E+08"
                                    ],
                                    "RCUT_NB": [
                                        "[angstrom] 1.0E+01"
                                    ]
                                }
                            ]
                        }
                    ],
                    "POISSON": [
                        {
                            "EWALD": [
                                {
                                    "ALPHA": [
                                        ".40"
                                    ],
                                    "EWALD_TYPE": [
                                        "NONE"
                                    ],
                                    "GMAX": [
                                        "80"
                                    ]
                                }
                            ],
                            "PERIODIC": [
                                "NONE"
                            ]
                        }
                    ]
                }
            ],
            "RESCALE_FORCES": [
                {
                    "MAX_FORCE": [
                        "10"
                    ]
                }
            ],
            "STRESS_TENSOR": [
                "NONE"
            ],
            "SUBSYS": [
                {
                    "CELL": [
                        {
                            "ABC": [
                                "[angstrom] 100 100 100"
                            ],
                            "ALPHA_BETA_GAMMA": [
                                "90 90 90"
                            ],
                            "PERIODIC": [
                                "NONE"
                            ]
                        }
                    ],
                    "TOPOLOGY": [
                        {
                            "CENTER_COORDINATES": [
                                {
                                    "CENTER_POINT": [
                                        "0 0 0"
                                    ]
                                }
                            ],
                            "CONN_FILE_FORMAT": [
                                "AMBER"
                            ],
                            "CONN_FILE_NAME": [
                                f"{TOP}"
                            ],
                            "COORD_FILE_FORMAT": [
                                "XYZ"
                            ],
                            "COORD_FILE_NAME": [
                                f"{XYZ}"
                            ]
                        }
                    ]
                }
            ]
        }
    ],
    "MOTION": [
        {
            "MD": [
                {
                    "VELOCITY_SOFTENING": [
                        {
                            "STEPS": [
                                f"{steps_vs}"
                            ],
                            "ALPHA": [
                                f"{alpha_vs}"
                            ],
                            "DELTA ": [
                                f"{delta_vs}"
                            ]
                        }
                    ],
                    "ENSEMBLE": [
                        "NVT"
                    ],
                    "STEPS": [
                        f"{md_steps}"
                    ],
                    "TEMPERATURE": [
                        f"{TEMPERATURE}"
                    ],
                    "THERMOSTAT": [
                        {
                            "CSVR": [
                                {
                                    "TIMECON": [
                                        "[fs] 0.5"
                                    ]
                                }
                            ],
                            "REGION": [
                                "GLOBAL"
                            ],
                            "TYPE": [
                                "CSVR"
                            ]
                        }
                    ],
                    "TIMESTEP": [
                        "[fs] 2.0"
                    ]
                }
            ],
            "PRINT": [
                {
                    "RESTART": [
                        {
                            "BACKUP_COPIES": [
                                "0"
                            ],
                            "EACH": [
                                {
                                    "MD": [
                                        "0"
                                    ]
                                }
                            ]
                        }
                    ],
                    "RESTART_HISTORY": [
                        {
                            "EACH": [
                                {
                                    "MD": [
                                        "0"
                                    ]
                                }
                            ]
                        }
                    ],
                    "TRAJECTORY": [
                        {
                            "EACH": [
                                {
                                    "MD": [
                                        "5"
                                    ]
                                }
                            ],
                            "FORMAT": [
                                "XYZ"
                            ]
                        }
                    ]
                }
            ]
        }
    ]
}

note = '''
- start from atropoisomer 1;
- total of 500 iterations (2 ns);
- all neighbours are exchanging;
- alphas are not constant: rise to a=0.05 if alpha <0.05;
-  NOT using external potential:
    C*((DIM)^2)^4, where C=10**-11; DIM = X,Y,Z
- potential affects only heavy atoms;
- no 'shaking' cells;
- default 1-4 interactions (1.0 electrostatic, 1.0 VdW);
- velocity softening: 1000 steps, alpha=0.01, default delta;
- T = 700 K;
- full contributions from torsions and improper torsions
'''

alpha_set = [0.0, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
# alpha_set = [0.0, 0.0001, 0.005, 0.01, 0.05, 0.1, 1.0]

# DIRNAME = 'ATROPOISOMERS-TEST-2-R-NL-HB-Q-RAB0-CN-POT-NOWALL'
# DIRNAME = 'Cu2L2-1'
DIRNAME = 'ATROPOISOMERS-VEL-SOFT-14'
XYZ_SHARED = 'structures.xyz'
EXCHANGE_STAT = 'exchanges.log'
ENERGY_STAT = 'energy.log'

NMAX = 500 # Total number of H-REMD cycles (500 for total MD time of 1 ns)
ALL_NEIGHBOURS = True # not tested yet
USE_LAST_STRUC = False # if lowest energy structure is in first half of trajectory, the last structure selected
# USE_WALL_POT = False
EXTERNAL_POT = False
USE_RMSD = False
CP2K_SHAKE_CELL = False

# 0. create Worlds
BETA = 1059.5151978
Nxyz = 0
try:
    os.mkdir(DIRNAME)
except:
    shutil.rmtree(DIRNAME, ignore_errors=True)
    os.mkdir(DIRNAME)
try:
    shutil.copyfile(XYZ_SHARED,os.path.join(DIRNAME,XYZ_SHARED))
    # get N of initial structures
    with open(XYZ_SHARED,'r') as fi:
        natoms = int(fi.readline())
        nlines = sum(1 for i in fi if i != '\n') + 1
        Nxyz = int(nlines/(2+natoms))
    print(f'{Nxyz} initial structures were read')
except:
    Nxyz = 0
    print('initial structures were not provided')
shutil.copyfile(XYZ,os.path.join(DIRNAME,XYZ))
# shutil.copyfile(PDB,os.path.join(DIRNAME,PDB))
shutil.copyfile(TOP,os.path.join(DIRNAME,TOP))
os.chdir(DIRNAME)
worlds = []
for alp in alpha_set:
    world = thread.World(
        BASE_DIR = os.path.abspath(os.getcwd()),
        alpha = alp,
        TOP = TOP,
        XYZ = XYZ,
        # XYZ0 = XYZ,
        XYZ_FOUND=XYZ_SHARED,
        ARGS = PARM, # __init__ takes deepcopy of dict
        USE_LAST_STRUC = USE_LAST_STRUC
    )
    worlds.append(world)
    if world.alpha < 0.001 and CP2K_SHAKE_CELL: # TODO: make it more clear
        world.shake = True
    if not EXTERNAL_POT:
        world.remove_ext_pot()

# start of loop
exchange_set = ['{:.3e}'.format(i) for i in alpha_set]
exchange_ord = {alp:exchange_set.index(alp) for alp in exchange_set}
with open(EXCHANGE_STAT,'w') as exc_stat:
    exc_stat.write('   '.join([str(exchange_ord[a]) for a in exchange_set]) + '\n')
with open(ENERGY_STAT,'w') as e_stat:
    e_stat.write('Type'.center(25) + ''.join([a.center(25) for a in exchange_set]) + '\n')
with open('comm.txt','w') as comm:
    comm.write(note)

world1 = worlds[-1]
#world1.remove_ext_pot()
N = 0
txtb = 0
t0 = time.time()
# if CP2K_SHAKE_CELL:
#     for w in worlds:
#         w.set_pot()

while N < NMAX:
    # 1. parallel launch of MDs
    procs = []
    main_queue = Queue()
    results = []
    for i,w in enumerate(worlds):
        if w.shake:
            w.shake_cell(N)
            # w.set_pot()
        calc = thread.MD_calc(w.XYZ,w.ARGS,w.CMD_LINE,w.INP,w.OUT,w.ERR,w.TRJ,w.alpha,Nxyz,w.XYZ_FOUND)
        proc = Process(target = calc.start_XTB,args=(main_queue,i))
        # proc = Process(target = thread2.start_XTB,args = (w.XYZ,w.ARGS,w.CMD_LINE,w.INP,w.OUT,w.ERR,w.TRJ,w.alpha,main_queue,i))
        procs.append(proc)
        proc.start()
    # loop to prevent deadlock, explained here:
    # stackoverflow.com/questions/31665328/
    stop_wait = False
    while not stop_wait:
        try:
            result = main_queue.get(False, 0.01)
            results.append(result)
        except queue.Empty:
            pass
        all_done = all([p.exitcode!=None for p in procs])
        stop_wait = all_done & main_queue.empty()
    results = sorted(results, key=lambda res: res['index'])
    for i,w in enumerate(worlds):
        w.set_results(results[i])
    txtb+=max([r['time']for r in results])
    # if all MDs failed, we exit
    if all([w.error for w in worlds]):
        print(f'All calculations are down, exiting (cycle {N})')
        break
    # 2. Metropolis–Hastings part
    for w in worlds:
        w.exchanged = False
    for i in range(len(worlds)-1):
        w1 = worlds[i]
        w2 = worlds[i+1]
        if (w1.error | w2.error):
            continue
        if w1.exchanged & ALL_NEIGHBOURS: # current implementation permits non-neighbors to exchange
            struc1 = w1.get_structure_after()
        else:
            struc1 = w1.get_structure_before()
        if w2.exchanged & ALL_NEIGHBOURS:
            struc2 = w2.get_structure_after()
        else:
            struc2 = w2.get_structure_before()
        alp1 = w1.alpha
        alp2 = w2.alpha
        e1 = struc1['Energy']
        e2 = struc2['Energy']
        rel_e = -e1*alp1 -e2*alp2 + e2*alp1 + e1*alp2 # delta = (e1-e2)*(alp2-alp1)
        if math.isnan(rel_e):
            continue
        delta = max(-400,min(400,rel_e * BETA)) # returns -400 if rel_e==nan
        p = math.exp(-delta)
        p0 = random.random()
        exchange = (p>=p0)
        print(f'exchange attempt: {alp1} and {alp2}')
        print(f'first index: {struc1["index"]}; second index: {struc2["index"]}')
        print(f'Relative energy and Delta: {rel_e} and {delta}')
        print(f'p and p0: {p} and {p0}')
        if USE_RMSD:
            struc1 = w1.get_structure_before()
            struc2 = w2.get_structure_before()
            val1 = struc1['RMSD']
            val2 = struc2['RMSD']
            exchange = (val1 < val2)
            print('Using RMSD criteria insted..')
            print(f'RMSD of {alp1}: {val1} | RMSD of {alp2}: {val2}')
        print(f'exchange: {exchange}')
        if exchange:
            w1.set_structure_after(struc2)
            w2.set_structure_after(struc1)
            w1.update_xyz()
            w2.update_xyz()
            w1.exchanged = True
            w2.exchanged = True
            # collect exchange statistics
            i1 = exchange_set.index('{:.3e}'.format(alp1))
            i2 = exchange_set.index('{:.3e}'.format(alp2))
            atmp = exchange_set[i1]
            exchange_set[i1] = exchange_set[i2]
            exchange_set[i2] = atmp
        else:
            w1.set_structure_after(struc1)
            w2.set_structure_after(struc2)
            w1.update_xyz()
            w2.update_xyz()
    with open(EXCHANGE_STAT,'a') as exc_stat:
        exc_stat.write('   '.join([str(exchange_ord[a]) for a in exchange_set]) + '\n')
    with open(ENERGY_STAT,'a') as e_stat:
        e_stat.write('ENERGY'.center(25) + ''.join(['{:.8e}'.format(w.struc_before_exchange['Energy']).center(25) for w in worlds]) + '\n')
        e_stat.write('VBIAS'.center(25) + ''.join(['{:.8e}'.format(w.struc_before_exchange['Vbias']).center(25) for w in worlds]) + '\n')
        if USE_RMSD:
            e_stat.write('RMSD'.center(25) + ''.join(['{:.8e}'.format(w.struc_before_exchange['RMSD']).center(25) for w in worlds]) + '\n')
    # 3. update shared XYZ
    world1.update_shared_xyz()
    Nxyz+=1
    
    # 4. prepare for next step
    for w in worlds:
        if w.error:
            #TODO: update .xyz with previos/from another world
            continue
        w.rename_trj_file()
    N+=1

    # 5. change order of exchange attempts
    # worlds.reverse() # TODO: better understand this feature

# postprocessing
postprocess.calc_rmsd(RMSD_REF,os.path.abspath(os.getcwd()),'rmsd.txt')
postprocess.calc_rmsd(RMSD_REF,os.path.abspath(os.getcwd()),'dihedral.txt',do_rmsd=False)

t = time.time()

if N==NMAX:
    print('Normal termination')
else:
    print(f'Error termination: all MDs failed at step {N}')
print('Total time: {:.3} h'.format((t-t0)/3600))
print('XTB time: {:.3} h'.format(txtb/3600))
print('Python overhead: {:.3} s'.format(t-t0-txtb))
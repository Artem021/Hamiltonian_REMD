import os, re, collections
import numpy as np

TEST = False

DEFAULT_PARM_CP2K = {
    'GLOBAL': {
        'SEED': '1642',
        'PROJECT': 'default', 
        'PRINT_LEVEL': 'MEDIUM', 
        'RUN_TYPE': 'MD'
    },
    'FORCE_EVAL': {
        'RESCALE_FORCES': {
            'MAX_FORCE': '10'
        }, 
        'STRESS_TENSOR': 'ANALYTICAL', 
        'METHOD': 'FIST', 
        'MM': {
            'FORCEFIELD': {
                'FORCE_SCALE': '1.0', 
                'EI_SCALE14': '1.0', 
                'ALL_VDW_SCALE': '1.0', 
                'ALL_EI_SCALE': '1.0', 
                'SHIFT_CUTOFF': 'F', 
                'PARMTYPE': 'AMBER', 
                'PARM_FILE_NAME': 'default.prmtop', 
                'SPLINE': {
                    'EMAX_SPLINE': '1.0E+08', 
                    'RCUT_NB': '[angstrom] 1.0E+01'
                }
            }, 
            'POISSON': {
                'EWALD': {
                    'EWALD_TYPE': 'NONE', 
                    'ALPHA': '.40', 
                    'GMAX': '80'
                }, 
                'PERIODIC': 'NONE'
            }
        }, 
        'SUBSYS': {
            'CELL': {
                'ABC': '[angstrom] 100 100 100', 
                'ALPHA_BETA_GAMMA': '90 90 90', 
                'PERIODIC': 'NONE'
                }, 
            'TOPOLOGY': {
                'CENTER_COORDINATES': {}, 
                'CONN_FILE_FORMAT': 'AMBER', 
                'CONN_FILE_NAME': 'default.prmtop', 
                'COORD_FILE_FORMAT': 'XYZ', 
                'COORD_FILE_NAME': 'default.XYZ'
            }
        }
    },
    'MOTION': {
        'MD': {
            'ENSEMBLE': 'NVT', 
            'TIMESTEP': '[fs] 2.0', 
            'STEPS': '1000', 
            'TEMPERATURE': '298', 
            'THERMOSTAT': {
                'REGION': 'GLOBAL', 
                'TYPE': 'CSVR', 
                'CSVR': {
                    'TIMECON': '[fs] 0.5'
                }
            }
        }, 
        'PRINT': {
            'RESTART': {
                'BACKUP_COPIES': '0', 
                'EACH': {
                    'MD': '0'
                }
            }, 
            'TRAJECTORY': {
                'FORMAT': 'XYZ', 
                'EACH': {
                    'MD': '5'
                }
            }, 
            'RESTART_HISTORY': {
                'EACH': {
                    'MD': '0'
                }
            }
        }
    }
}

DEFAULT_PARM = {
    '$md':{'restart':'false','hmass':4,'time':10.00,'temp':300.00,'step':2.00,'shake':0,'dump':10.00},
    '$metadyn':{'save':5,'kpush':0.042000,'alp':1.300000,'alpha':0.0},
    '$constrain':{'force constant':0.25,'all bonds':'T','reference':'reference.xyz'},
    '$wall':{'potential':'logfermi','sphere':'auto,all'}
}

SPECIAL_SIGNS = {'elements':':','distance':':','sphere':':'}

#os.chdir('cp2k')

def read_cp2k(fname,fi=None):
    sect = {}
    if fi == None:
        fi = open(fname,'r')
    while True:
        line = fi.readline()
        if line == '\n':
            continue
        if 'END' in line or line == '':
            return sect
        line = line.strip()
        if '&' in line:
            sctname = line[1:]
            subsect = read_cp2k(fname,fi)
            if sect.get(sctname) != None:
                sect[sctname].append(subsect)
            else:
                sect[sctname] = [subsect]
            continue
        key,val = line.split(maxsplit=1)
        if sect.get(key) != None:
            sect[key].append(val)
        else:
            sect[key] = [val]

def write_cp2k(sect_dict):
    lines = []
    def write_rec(sect_dict,l):
        for f in sect_dict:
            for i in sect_dict[f]:
                if type(i) == dict:
                    lines.append(' '*l+'&'+f)
                    l+=2
                    write_rec(i,l)
                    l-=2
                    lines.append(' '*l+'&END '+f)
                else:
                    if f == '':
                        lines.append(''*l+i.strip())
                    else:
                        lines.append(' '*l+f+' '+i)
        return '\n'.join(lines)+'\n'
    return(write_rec(sect_dict,1))

# import json

# sections = read_cp2k('cp2k.inp')
# print(sections)

# print(sections['FORCE_EVAL'][0]['MM'][0]['FORCEFIELD'][0]['FORCE_SCALE'][0])
# print(json.dumps(sections,sort_keys=True, indent=4))

# exit()

# with open('cp2k2.inp','w') as fo:
#     txt = write_cp2k(sections)
#     fo.write(txt)

def read_section(template_file_name,sect_name):
    with open(template_file_name,'r') as tm:
        sect_name_lines = tm.readlines()
    for i,line in enumerate(sect_name_lines):
        if sect_name in line:
            start = i
            if not '&' in sect_name:
                assert '&' not in line
            break
    def read_recursively(s):
        x = {}
        lock = 0
        skip = False
        for i in range(s,len(sect_name_lines)):
            if '&' in sect_name_lines[i]:
                if 'END' in sect_name_lines[i]:
                    lock -=1
                    if lock <= 0:
                        skip = True
                else:
                    if lock ==0 and skip ==False:
                        res = re.search('([A-Z|0-9|_]+)',sect_name_lines[i])
                        sub_sect_name = res.group(1)
                        x[sub_sect_name] = {}
                    elif lock ==1 and skip ==False:
                        for f in read_recursively(i): x[sub_sect_name][f] = read_recursively(i)[f]
                    lock +=1
            else:
                if lock <=1 and skip ==False:
                    pattern = re.compile('^\s+([A-Z|_|0-9]+) (.*)')
                    res = re.search(pattern,sect_name_lines[i])
                    if res == None:
                        if '' in x[sub_sect_name]:
                            x[sub_sect_name]['']+=sect_name_lines[i]
                        else:
                            x[sub_sect_name][''] = sect_name_lines[i]
                    else:
                        x[sub_sect_name][res.group(1)] = res.group(2)
        return x
    if '&' in line:
        return read_recursively(start)
    else:
        pattern = re.compile('({}) (.*)'.format(sect_name))
        res = re.search(pattern,line)
        return({res.group(1):res.group(2)})

def write_section(sect_dict):
    lines = []
    def write_recursively(sect_dict,l):
        for f in sect_dict:
            if type(sect_dict[f]) == dict:
                lines.append(' '*l+'&'+f)
                l+=2
                write_recursively(sect_dict[f],l)
                l-=2
                lines.append(' '*l+'&END '+f)
            else:
                if f == '':
                    lines.append(''*l+sect_dict[f].strip())
                else:
                    lines.append(' '*l+f+' '+sect_dict[f])
        return '\n'.join(lines)+'\n'
    return(write_recursively(sect_dict,1))

def write_input_cp2k(parm,fname='cp2k.inp'):
    with open(fname,'w') as fo:
        fo.write(write_cp2k(parm))
        # fo.write(write_section(parm))




# os.chdir('cp2k')
# glob = read_section('cp2k.inp','&GLOBAL')
# mot = read_section('cp2k.inp','&MOTION')
# fe = read_section('cp2k.inp','&FORCE_EVAL')
# print(fe)
# write_input_cp2k(DEFAULT_PARM_CP2K,fname = 'cp2k0.inp')








def reverse_readline(filename, buf_size=8192):
    """A generator that returns the lines of a file in reverse order"""
    with open(filename) as fh:
        segment = None
        offset = 0
        fh.seek(0, os.SEEK_END)
        file_size = remaining_size = fh.tell()
        while remaining_size > 0:
            offset = min(file_size, offset + buf_size)
            fh.seek(file_size - offset)
            buffer = fh.read(min(remaining_size, buf_size))
            remaining_size -= buf_size
            lines = buffer.split('\n')
            # The first line of the buffer is probably not a complete line so
            # we'll save it and append it to the last line of the next buffer
            # we read
            if segment is not None:
                # If the previous chunk starts right from the beginning of line
                # do not concat the segment to the last line of new chunk.
                # Instead, yield the segment first
                if buffer[-1] != '\n':
                    lines[-1] += segment
                else:
                    yield segment
            segment = lines[0]
            for index in range(len(lines) - 1, 0, -1):
                if lines[index]:
                    yield lines[index]
        # Don't yield None if the file was empty
        if segment is not None:
            yield segment

def update_parm(new,base=DEFAULT_PARM):
    sections = new.keys() & base.keys()
    assert len(sections.difference(new.keys()))==0, f"unknown section(s) in: {new}"
    for section in sections:
        base[section].update(new[section])
    return base
    

def write_input(parm,fname='xtb.inp'):
    with open(fname,'w') as fo:
        for section in parm:
            fo.write(section+'\n')
            subsect = parm[section]
            for key in subsect:
                sign = SPECIAL_SIGNS.get(key,'=')
                fo.write(f'  {key}{sign}{subsect[key]}\n')
        fo.write('$end')

def parse_input(fname='xtb.inp'):
    parm = {}
    keys = {}
    section = None
    with open(fname,'r') as fi:
        lines = fi.readlines()
    delim = set(SPECIAL_SIGNS.values())
    delim.add('=')
    pat = f'[{"".join(delim)}]'
    # pat = re.compile(f'\s+([a-z]+)\s*[{",".join(delim)}]\s*(.*)')
    for line in lines:
        if '$' in line:
            if section!=None:
                parm[section] = keys.copy()
            section = line.strip()
            keys.clear()
        else:
            kw,val = re.split(pat,line.strip())
            keys[kw] = val
    return parm

def get_frame_xyz(fname,index=-1,frames=None,as_np=False,nmax=1000):
    xyz = []
    na = 0
    if index == -1:
        lines = iter(reverse_readline(fname))
        while True:
            line = next(lines).strip()
            # line = line.strip()
            if 'xtb' in line or 'time' in line:
                comment = line
                break
            if na > nmax:
                raise RuntimeError(f'incorrect .xyz file: {fname}, check comment line \
                    (must contain "xtb") and N atoms (if N > 1k, add "nmax = N" to call)')
            xyz.append(line)
            na+=1
        na0 = next(lines).strip()
        assert int(na0)==na, f'incorrect N atoms: {na} was read, expected {na0}'
        if as_np:
            xyz.reverse()
            xyz_np = np.array([l.split()[1:] for l in xyz],dtype=float)
            return xyz_np
        else:
            xyz.append(comment)
            xyz.append(na0)
            xyz.reverse()
            return '\n'.join(xyz)
    else:
        assert index>=0
        with open(fname,'r') as fi:
            nframes = 0
            eof = False
            line1 = fi.readline().strip()
            nat = int(line1)
            buf = collections.deque([line1],nat+3)
            while True:
                line = fi.readline().strip()
                if 'xtb' in line or 'time' in line:
                    nframes+=1
                if nframes > index+1:
                    break
                if line == '':
                    eof = True
                    break
                buf.append(line.strip())
            assert nframes > 0, f'incorrect .xyz file: {fname}, check comment line (must contain "xtb")'
            if eof:
                buf.popleft()
            else:
                buf.pop()
            if as_np:
                buf.popleft() # - natoms
                buf.popleft() # - comment
                xyz_np = np.array([l.split()[1:] for l in buf],dtype=float)
                return xyz_np
            else:
                return '\n'.join(buf)

# xyz = get_frame_xyz('example/xtb.trj',-1)
# with open('test-m.xyz','w') as fo:
#     fo.write(xyz)
#     fo.write('\n')
#     fo.write(xyz)

# write_input(DEFAULT_PARM,'xtb.inp')
# d = parse_input()
# write_input(d,'xtb2.inp')


if TEST:
    import time
    #xyz (eof)
    start = time.time()
    xyz = get_frame_xyz('example/xtb.trj',-1)
    end = time.time()
    print(f'from EOF: {end-start} seconds')
                
    #xyz (bof)
    start = time.time()
    xyz = get_frame_xyz('example/xtb.trj',999)
    end = time.time()
    print(f'from BOF: {end-start} seconds')

    #xyz numpy
    xyz_np = get_frame_xyz('example/xtb.trj',as_np=True)
    print(xyz_np)

    # write/parse
    write_input(DEFAULT_PARM,'xtb.inp')
    d = parse_input()
    write_input(d,'xtb2.inp')
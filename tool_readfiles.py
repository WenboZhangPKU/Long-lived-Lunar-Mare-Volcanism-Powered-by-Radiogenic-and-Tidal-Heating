import numpy as NP
import os
import subprocess

#----read in radius----#
def read_radius(in_dict,route,casename):
    nprocz = get_variable(in_dict,'nprocz','int')
    noz = get_variable(in_dict,'nodez','int')
    lnoz=int((noz-1)/nprocz)+1
    data = read_horig_output(route,casename,0,nprocz,lnoz,1)
    return data[:,0]

#----read in horizontal average result---#
def read_horig_output(Croute,Ccasename,Cstep,Cnprocz,Clnoz,cols):
    noz = (Clnoz-1)*Cnprocz+1
    Rdata = NP.zeros((noz,cols))
    for pz in range(Cnprocz):
        filename = "%s/%s.horiz_avg.%d.%d"%(Croute,Ccasename,pz,Cstep)
        fin = open(filename,'r')
        offset = (Clnoz-1)*pz
        nl = 0
        for line in fin:
            temp = line.split(' ')
            for nc in range(cols):
                Rdata[nl+offset][nc] = float(temp[nc])
            nl = nl+1
        fin.close()
    return Rdata

#----read in time----#
def read_time(Croute,Ccasename):
    filename = "%s/%s.time"%(Croute,Ccasename)
    fin = open(filename,'r')
    step = list()
    time = list()
    for line in fin:
        temp = line.split(' ')
        step.append(int(temp[0]))
        time.append(float(temp[1]))
    #----drop duplicate data----#
    fstep = step[len(step)-1]
    Rdata = NP.zeros(fstep+1)
    temp = fstep+1
    for n in range(len(step)-1,-1,-1):
        n_s = step[n]
        if n_s<temp:
            Rdata[n_s] = time[n]
            temp = n_s
    fin.close();
    return Rdata
#-----read in step and time-----#
#Rdata are completely filled, that is, each row of Rdata is filled#
def read_time2(Croute,Ccasename):
    filename = "%s/%s.time"%(Croute,Ccasename)
    fin = open(filename,'r')
    step = list()
    time = list()
    for line in fin:
        temp = line.split(' ')
        step.append(int(temp[0]))
        time.append(float(temp[1]))
    #----drop duplicate data----#
    fstep = step[len(step)-1]
    Rdata = NP.zeros((fstep+1,2))
    temp = fstep+1
    for n in range(len(step)-1,-1,-1):
        n_s = step[n]
        if n_s<temp:
            Rdata[n_s,0] = step[n]
            Rdata[n_s,1] = time[n]
            temp = n_s
    fin.close();
    return Rdata

#-----read in heatflux-----#
#Rdata usually has spare rows#
def read_heatflux(Croute,Ccasename):
    filename = "%s/%s.qb.dat"%(Croute,Ccasename)
    fin = open(filename,'r')
    time = list()
    heatflux =list()
    sqrtvdotv =list()
    Tcmb =list()
    rin =list()
    for line in fin:
        temp = line.split('  ')
        time.append(float(temp[1]))
        heatflux.append(float(temp[2]))
        sqrtvdotv.append(float(temp[3]))
        Tcmb.append(float(temp[4]))
        rin.append(float(temp[5]))
    #-----drop redundant data-----#
    Rdata = NP.zeros((len(time),5))
    ftime =time[len(time)-1]+1e-9
    n_s =len(time)-1
    for i in range(len(time)-1,-1,-1):
        if time[i]<ftime:
            Rdata[n_s,0] = time[i]
            Rdata[n_s,1] = heatflux[i]
            Rdata[n_s,2] = sqrtvdotv[i]
            Rdata[n_s,3] = Tcmb[i]
            Rdata[n_s,4] = rin[i]
            n_s =n_s-1
            ftime =time[i]+1e-9
    fin.close();
    return Rdata

#----read in q_file----#
def read_q_file(filename,time_array,column,colt=0):
    tiny = 1e-4
    ttiny = 1e-9
    #----read data from qfile----#
    fin = open(filename,'r')
    data = read_data(filename)
    lines = len(data)//column
    #----eliminate overlapping data----#
    step = list()
    ilist = list()
    time_after = 1.1*data[column*(lines-1)+colt]
    for i in range(lines-1,-1,-1):
        time = data[column*i+colt]
        if time<time_after:
            ilist.append(i)
            time_after = time
            if time<ttiny:
                step.append(0)
                continue
            for n in range(len(time_array)):
                if abs(time-time_array[n])/time<tiny:
                    step.append(n)
                    break
    #----generate data matrix for return----#
    step = step[::-1]
    Rdata = NP.zeros((len(ilist),column))
    for n in range(len(ilist)-1,-1,-1):
        i = ilist[n]
        for j in range(column):
            Rdata[n,j] = data[column*i+j]
    return Rdata,step

#----read input_file----#
def read_input(Cfilename):
    Rdict = dict()
    fin = open(Cfilename,'r')
    for line in fin:
        #----elimate space at head----#
        for i in range(len(line)):
            if line[i]!=' ' and line[i]!='\t':
                break
        line=line[i:len(line)]
        #----skip some lines----#
        if (line is '\n') or (line[0] is '#'):
            continue
        #----elimate content after '#'----#
        for i in range(len(line)):
            if line[i] is '#':
                break
        line=line[0:i]
        temp = line.split('=')
        if temp[1][0] is '"':
            temp[1] = temp[1][1:len(temp[1])-1]
        #print(line) #debug
        Rdict[temp[0]]=temp[1]
    fin.close()
    return Rdict

#----get value of a variable from an existing dictionary----#
def get_variable(Cdict,Cname,Ctype):
    try:
        temp=Cdict[Cname]
    except NameError:
        print('no name %s in dictionary\n',Cname)
        os.exit()
    #print(temp) #debug
    if Ctype is 'string':
        Rvalue=temp
    elif Ctype is 'int':
        Rvalue=int(temp)
    elif Ctype is 'float':
        Rvalue=float(temp)
    elif Ctype is 'int_list':
        Rvalue=temp.split(',')
        for i in range(len(Rvalue)):
            Rvalue[i] = int(Rvalue[i])
    elif Ctype is 'float_list':
        Rvalue=temp.split(',')
        for i in range(len(Rvalue)):
            Rvalue[i] = float(Rvalue[i])
    return Rvalue

def read_field(route,case_name,nprocz,noz,ll_max,step):
    header = 5
    inter = 2
    columns = 6
    pcolumns = 2
    Dname = ('T')
    f_dict = {'T':[2,3],'Vr':[4,5]} #columns, starts from 0
    lnoz=int((noz-1)/nprocz)+1
    nn_max = int((ll_max+1)*(ll_max+2)/2)
    Rdata = NP.zeros((noz,ll_max+1)) #
    for procz in range(nprocz):
        f_file = "%s/%s.field.%d.%d"%(route,case_name,procz,step)
        data = read_data(f_file)
        for i in range(lnoz):
            nz = i+(lnoz-1)*procz
            Atemp = NP.zeros(ll_max+1)
            for ll in range(ll_max+1):
                for mm in range(ll+1):
                    pp = int(ll*(ll+1)/2+mm)
                    for name in Dname:
                        idx0 = header+(i*nn_max+pp)*columns+inter*i+f_dict[name][0]
                        idx1 = header+(i*nn_max+pp)*columns+inter*i+f_dict[name][1]
                        Atemp[ll] = Atemp[ll]+data[idx0]**2.0+data[idx1]**2.0
                Atemp[ll] = Atemp[ll]**0.5 ## It refers to the total power of the field at spherical harmonic degree l, 
                                           ##   which we will call the power per degree, from Wieczorek & Meschede, 2018
            Rdata[nz,:] = Atemp
    return Rdata

def write_surf_runfile(in_dict,route,case_name,step,vtype='total'):
    #----patameters----#
    ofile = './cc/runfile'
    caps = 12
    nox = get_variable(in_dict,'nodex','int')
    noy = get_variable(in_dict,'nodey','int')
    noz = get_variable(in_dict,'nodez','int')
    nprocx = get_variable(in_dict,'nprocx','int')
    nprocy = get_variable(in_dict,'nprocy','int')
    nprocz = get_variable(in_dict,'nprocz','int')
    inputf0="%s/%s"%(route,case_name) #
    outputf0="./cc/%s"%(case_name) # output result in the same as input
    cpu_xy = caps*nprocx*nprocz
    cpu_z = nprocz
    file_type=1
    get_deltaT=3
    comp_temp=1
    get_slab_center=40
    slab_t=-0.1
    special_value=1.0
    #----output to runfile----#
    with open(ofile,'w') as fout:
        fout.write('%s\n'%(inputf0))
        fout.write('%s\n'%(inputf0))
        fout.write('%s\n'%(outputf0))
        fout.write('%d %d %d\n'%(nox,noy,noz))
        fout.write('%d %d\n'%(cpu_xy,cpu_z))
        fout.write('%d %d %d 0 %d\n'%(file_type,get_deltaT,step,comp_temp))
        fout.write('%d %.4f 0.0 %.4f\n'%(get_slab_center,slab_t,special_value))
        fout.write('-1')


#----read from a pure data file----#
def read_data(filename):
    data=[]
    with open(filename) as f:
        for line in f:
            line.strip()
            for part in line.split():
                if part is not ' ':
                    data.append(float(part))
    return data
#----read data and return a matrix----#
def read_data_1(filename):
    with open(filename,'r') as fin:
        line=fin.readline()
        s_list=line.split(' ')
    col=0
    for s in s_list:
        if s == ' ' or s == '\n' or s == '':
            pass
        else:
            col=col+1
    in_data=read_data(filename)
    row=len(in_data)//col
    data=NP.array(in_data)
    data.resize((row,col))
    return data
#----get steps by finding horiz file----#

def get_steps(Croute,Ccase_name,Cmax_step=200001):
    Rtuple = ()
    for step in range(Cmax_step):
        filename = "%s/%s.horiz_avg.0.%d"%(Croute,Ccase_name,step)
        if os.path.isfile(filename):
            Rtuple=Rtuple+(step,)
    return Rtuple

#----write to file when file or variable doesn't exit----#
def write_to_file1(pp_file,name,value,vtype='int',vvtype='simple'):
    if os.path.isfile(pp_file):
        pp_dict = read_input(pp_file)
        try:
            pp_dict[name]
        except KeyError:
            write = 1
        else:
            write = 0
    else:
        write = 1
    if write:
        write_to_file(pp_file,name,value,vtype,vvtype)

def overwrite_to_file(pp_file,name,value,vtype='int',vvtype='simple'):
    temp_file = '%s_temp'%(pp_file)
    if os.path.isfile(pp_file):
        pp_dict = read_input(pp_file)
    else:
        pp_dict = []
    for key in pp_dict:
        if key == name:
            write_to_file(temp_file,name,value,vtype,vvtype)
        else:
            try:
                fp = open(temp_file,'a')
            except IOError:
                print('cannot open %s for writing',pp_file)
            else:
                fp.write('%s=%s\n'%(key,pp_dict[key]))
                fp.close()
    os.rename(temp_file,pp_file)

#----write to post-process file----#
def write_to_file(pp_file,name,value,vtype='int',vvtype='simple',vvvtype='value'):
    with open(pp_file,'a') as fout:
        if vtype is 'int':
            s='%d'
        elif vtype is 'float':
            s='%.4e'
        elif vtype is 'str':
            s='%s'
        if vvtype is 'array':
            string = ''
            for i in range(len(value)):
                string = string+s%(value[i])
                if i is not len(value)-1:
                    string = string+','
        else:
            string = s%value
        if vvvtype is 'quote':
            string = '\''+string+'\''
        fout.write("%s="%(name)+string+'\n')

#----run bash file----#
def run_bash(bashcommand,dcwd=None,shll=False):
    #----output output to a specific file----#
    for command in bashcommand:
        if shll is True:
            process = subprocess.Popen(command, stdout=subprocess.PIPE, cwd=dcwd,shell=True)
        else:
            process = subprocess.Popen(command.split(), stdout=subprocess.PIPE,cwd=dcwd)
        output, error = process.communicate()
    return

#----read data from vtk files, return as data and coordinate as numpy arrays----#
def read_vtk_file(filename):
    def _data_type(line):
        d_type={}
        d_list=line.split(' ')
        name=None
        comp=1
        for string in d_list:
            index=string.find('=')
            if index > -1:
                head=string[:index]
                value=string[index+2:len(string)-1]
                if head == 'Name':
                    name=value
                elif head == 'NumberOfComponents':
                    comp=int(value)
        return name,comp
    d_dict={}
    with open(filename,'r') as fin:
        for i in range(5):
            line=fin.readline()
        while 1:
            line=fin.readline()
            name,comp=_data_type(line)
            if name==None:
                break
            data=[]
            while 1:
                line=fin.readline()
                if line.find('<') > -1:
                    break
                d_list=line.split(' ')
                length=len(d_list)
                for j in range(length):
                    data.append(float(d_list[j]))
            data=NP.array(data)
            row=len(data)//comp
            data.resize((row,comp))
            d_dict[name]=data
        #read coordinate#
        data=[]
        comp=3
        for i in range(4):
            line=fin.readline()
        while 1:
            line=fin.readline()
            if line.find('<') > -1:
                break
            d_list=line.split(' ')
            length=len(d_list)
            for j in range(length):
                data.append(float(d_list[j]))
        data=NP.array(data)
        row=len(data)//comp
        data.resize((row,comp))
        coord=data
    return d_dict,coord
#----write vtk file from data and coord, should be numpy arrays----#
def write_vtk_file(filename,o_data,o_coord,step,proc,nx,ny,nz):
    with open(filename,'w') as fin:
        #vtk head
        fin.write('<?xml version=\"1.0\"?>\n')
        fin.write('<VTKFile type=\"StructuredGrid\" version=\"0.1\" \
compressor=\"vtkZLibDataCompressor\" byte_order=\"LittleEndian\">\n')
        fin.write('  <StructuredGrid WholeExtent="%d %d %d %d %d %d">\n'\
                %(nz[0],nz[1],nx[0],nx[1],ny[0],ny[1]))
        fin.write('    <Piece Extent="%d %d %d %d %d %d">\n'\
                %(nz[0],nz[1],nx[0],nx[1],ny[0],ny[1]))
        fin.write('      <PointData Scalars=\"temperature\" Vectors=\"velocity\">\n')
        #vtk data
        for key in o_data:
            data=o_data[key]
            if data.shape[1]>1:
                comp=data.shape[1]
                fin.write('        <DataArray type=\"Float32\" Name=\"%s\"\
 NumberOfComponents=\"%d\" format=\"ascii\">\n'%(key,comp))
            else:
                fin.write('        <DataArray type=\"Float32\" Name=\"%s\" format=\"ascii\">\n'%(key))
            for i in range(data.shape[0]):
                for j in range(data.shape[1]):
                    if j == data.shape[1]-1:
                        fin.write('%.4e\n'%(data[i,j]))
                    else:
                        fin.write('%.4e '%(data[i,j]))
            fin.write('        </DataArray>\n')
        fin.write('      </PointData>\n')
        fin.write('      <CellData>\n')
        fin.write('      </CellData>\n')
        fin.write('      <Points>\n')
        #vtk coordinate
        fin.write('        <DataArray type=\"Float32\" Name=\"coordinate\"\
 NumberOfComponents=\"3\" format=\"ascii\">\n')
        for i in range(o_coord.shape[0]):
            for j in range(o_coord.shape[1]):
                if j == o_coord.shape[1]-1:
                    fin.write('%.4e\n'%(o_coord[i,j]))
                else:
                    fin.write('%.4e '%(o_coord[i,j]))
        fin.write('        </DataArray>\n')
        fin.write('     </Points>\n')
        fin.write('    </Piece>\n')
        fin.write('  </StructuredGrid>\n')
        fin.write('</VTKFile>\n')

#----------------------------------------------------------------#
#find the very time step we wanna plot#
#only useful in continuous from Step 0 to Stepxxx#
#regardless of changed storage_spacing#
#----------------------------------------------------------------#
def find_timestep(route,case_name,time_array,t_real,thermdiff,radius):
    #route = "/home/p300/newdisk/back/input/D150v1e-1V1e21H100_20201224"
    #case_name = 'a'
    #time_array =readf.read_time(route,case_name)
    #t_real =4.7e8#4.42e9 starting time
    #t_real =6.2e8
    t_nondimensionalized =thermdiff*t_real*365*24*60*60/radius**2
    t_position =0
    diff_temp =t_nondimensionalized
    for  i in  range(1,len(time_array)):
        if abs(time_array[i] -t_nondimensionalized)< abs(diff_temp):
            t_position =i
            diff_temp =time_array[i] -t_nondimensionalized
    lowerbound=int(NP.floor(t_position/100)*100)
    time_lower=time_array[lowerbound]*radius**2/thermdiff/(365*24*60*60)/1e6
    upperbound=lowerbound +100
    time_upper=time_array[upperbound]*radius**2/thermdiff/(365*24*60*60)/1e6
    print("Timestep=%d, Timestep%d=%.5g Myr, Timestep%d=%.5g Myr \n" %(t_position,lowerbound,time_lower,upperbound,time_upper))

#only useful for nprocx=nprocy=nprocz=2#
def count_upwellingpoint(solution_cycle_init,nproc,nno,noz,route,IBCthickness,
                        lnox=33, lnoy=33):
    vtk_path =route
    point_upwelling =0
    point_boundary =0 ##all the boundary points (vertex+ other boundary points)
    point_vertex3 =0 ## vertex of contribution =3
    point_vertex4 =0 ## vertex of contribution =4
    tiny =1e-8
    lmesh = lnox * lnoy
    for i in range(int(nproc/2)):
        n =2*i+1
        vtk_up ='a.proc'+str(n)+'.'+str(solution_cycle_init)+'.vts' ##450-km-depth, nr=36 for 150km; 38 for 100km
        vtk_up =os.path.join(vtk_path,vtk_up)
        #read in composition from upper and lower blocks#
        data_up, coord_up =read_vtk_file(vtk_up)
        comp1_up =data_up['composition1']
        velocity_up=data_up['velocity']
        comp1_threshold =0.4 ##
        if IBCthickness ==150:
            nr_450kmdepth=36
            comp1_threshold =comp1_threshold
        elif IBCthickness ==100:
            nr_450kmdepth=38
            comp1_threshold =comp1_threshold*42.57/61.95
        elif IBCthickness ==75:
            nr_450kmdepth=40
            comp1_threshold =comp1_threshold*42.57/81.37
        elif IBCthickness ==50:
            nr_450kmdepth=40
            comp1_threshold =comp1_threshold*42.57/120.2
        else:
            print("Thickness doesn't match!")
        for j in range(lmesh): ##33*33 points on a 450-km-depth small cap
            if (comp1_up[j*33+nr_450kmdepth-noz]>comp1_threshold): ##try
                velocity_radial=0
                for k in range(3):
                    velocity_radial +=coord_up[j*33+nr_450kmdepth-noz][k]*velocity_up[j*33+nr_450kmdepth-noz][k]
                if velocity_radial>0: ## separated from Line 500 to reduce computations
                    if j>=0 and j<=32: ## the first edge, anticlockwise
                        if j==0: ## the first vertex
                            if n==9 or n==33 or n==57 or n==81: ## the 9th, 33rd, 57th and 81st cap
                                point_vertex3 +=1 ## contribution=3
                            else:
                                point_vertex4 +=1
                        elif j==32: ## the second vertex
                            if n==3 or n==27 or n==51 or n==75 or n==19 or n==43 or n==67 or n==91:
                                point_vertex3 +=1 ## contribution=3
                            else:
                                point_vertex4 +=1
                        else:
                            point_boundary +=1 ## other edeg points
                    elif j>=33*32 and j<=33*33-1: ## the third edge, anticlockwise
                        if j==33*32: ## the forth vertex
                            if n==5 or n==29 or n==53 or n==77 or n==21 or n==45 or n==69 or n==93:
                                point_vertex3 +=1
                            else:
                                point_vertex4 +=1
                        elif j==33*33-1: ## the third vertex
                            if n==15 or n==39 or n==63 or n==87:
                                point_vertex3 +=1
                            else:
                                point_vertex4 +=1
                        else:
                            point_boundary +=1 ## other edge points
                    elif j%33 ==0: ## the forth edge
                        point_boundary +=1
                    elif j%33 ==32: ## the second edge
                        point_boundary +=1
                    else:
                        point_upwelling +=1 ## internal points
    point_vertex3 =point_vertex3/3
    point_vertex4 =point_vertex4/4
    point_boundary =point_boundary/2
    point_upwelling =point_upwelling +point_vertex3 +point_vertex4 +point_boundary
    ratio =(point_upwelling)/(46128+8+42+2976)
    print("\n step=%d" %solution_cycle_init)
    print("\n ratio=%.5e" %ratio)
    return ratio

def binary_search(ratio_threshold,head,tail,nproc,nno,noz,route,IBCthickness):
    mid =int((head+tail)/2)
    #tiny =1e-8
    ratio =count_upwellingpoint(mid*100, nproc, nno, noz, route, IBCthickness)
    if abs(ratio-ratio_threshold) <1.0e-4:
        return mid,ratio
    elif tail-head <=1:
        return tail,ratio
    elif ratio <ratio_threshold:
        binary_search(ratio_threshold, mid, tail, nproc, nno, noz, route, IBCthickness)
    elif ratio >ratio_threshold:
        binary_search(ratio_threshold, head, mid, nproc, nno, noz, route, IBCthickness)

def count_point(solution_cycle_init,nproc,nno,noz,route,IBCthickness):
    vtk_path =route
    point_total =0
    point_upwelling =0
    point_boundary =0 ##total
    point_boundary1 =0 ##vertex
    #point_duplicate =0 ##total
    #point_duplicate1 =0 ##vertex
    for i in range(int(nproc/2)):
        vtk_up ='a.proc'+str(2*i+1)+'.'+str(solution_cycle_init)+'.vts' ##450-km-depth, nr=36 for 150km; 38 for 100km
        vtk_up =os.path.join(vtk_path,vtk_up)
        #read in composition from upper and lower blocks#
        data_up, coord_up =read_vtk_file(vtk_up)
        comp1_up =data_up['composition1']
        velocity_up=data_up['velocity']

        if IBCthickness ==150:
            nr_450kmdepth=36
        elif IBCthickness ==100:
            nr_450kmdepth=38
        elif IBCthickness ==50:
            nr_450kmdepth ==40
        else:
            print("Thickness doesn't match!")
        for j in range(33*33):
            point_total +=1
            if j>=0 and j<=32:
                point_boundary +=1
                if j==0 or j==32:
                    point_boundary1 +=1
            elif j>=33*32 and j<=33*33-1:
                point_boundary +=1
                if j==33*32 or j==33*33-1:
                    point_boundary1 +=1
            elif j%33 ==0:
                point_boundary +=1
            elif j%33 ==32:
                point_boundary +=1
    #ratio =(point_upwelling-(point_duplicate-point_duplicate1)/2-point_duplicate1*0.75)/(52272-3552)
    print("\n total contribution of surface points =%d" %point_total)
    #print("\n ratio=%.5e" %ratio)
    #print("\n step_upwelling =%d" %(solution_cycle_init))
    print("\n contribution of boundary points =%d" %point_boundary)
    print("\n contribution of vertex points =%d" %point_boundary1)
    #return ratio

#----write vtm file ----#
def write_vtm_file(filename,step,nproc):
    with open(filename,'w') as fin:
        #vtm head
        fin.write('<?xml version=\"1.0\"?>\n')
        fin.write('<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" \
compressor=\"vtkZLibDataCompressor\" byte_order=\"LittleEndian\">\n')
        fin.write('  <vtkMultiBlockDataSet>\n')
        #vtm file
        for i in range(nproc):
            fin.write('    <DataSet index=\"%d\" file=\"newa.proc%d.%d.vts\"/>\n' %(i,i,step))
        #vm tail
        fin.write('  </vtkMultiBlockDataSet>\n')
        fin.write('</VTKFile>\n')

def read_vtk_file_extent(filename):
    with open(filename,'r') as fin:
        for i in range(2):
            line=fin.readline()
        line=fin.readline()
        extent=line.split('\"')[1]
        extent =extent.split( )
        nz =[int(extent[0]),int(extent[1])]
        nx =[int(extent[2]),int(extent[3])]
        ny =[int(extent[4]),int(extent[5])]
    return nz,nx,ny

#-------------read profiles of initial condition-------------#
#such as nondimensionalized temperature profile, pressure profile#
#solidus and liquidus#
def read_ic_file(filename):
    with open(filename,'r') as profile:
        n=0
        radius =[]
        ndvalue =[]
        for row in profile.readlines():
            row =row.split(' ')
            radius.append(float(row[0]))
            ndvalue.append(float(row[1]))
            n +=1
        return radius,ndvalue

#----read in hemispherically horizontal average result---#
def read_hemishorig_output(Croute,Ccasename,Cstep,Cnprocz,Clnoz,cols):
    noz = (Clnoz-1)*Cnprocz+1
    Rdata = NP.zeros((noz,cols))
    for pz in range(Cnprocz):
        filename = "%s/%s.hemishoriz_avg.%d.%d"%(Croute,Ccasename,pz,Cstep)
        fin = open(filename,'r')
        offset = (Clnoz-1)*pz
        nl = 0
        for line in fin:
            temp = line.split(' ')
            for nc in range(cols):
                Rdata[nl+offset][nc] = float(temp[nc])
            nl = nl+1
        fin.close()
    return Rdata

#---get paths of input files---#
def get_input_file(route,input_route):
    name=route.split('/')[-1]
    if name=='':
        name=route.split('/')[-2]
    serialnumber=name.split('_')[-1]
    print(serialnumber)
    input_file = input_route+"in_65moon_"+serialnumber
    try:
        thickness=int(name[1:4])
    except ValueError:
        thickness=int(name[1:3])
    return input_file, serialnumber, thickness

 #----creat dirname if not exit----#
def assign_output_dir(*args):
     i = 0
     for mem in args:
         if i == 0:
             dirname = mem
         else:
             dirname = os.path.join(dirname,mem)
         i = i+1
     if os.path.isdir(dirname) is False:
         os.mkdir(dirname)
     return dirname

#----------------------------------------------#
#get one type of the field data at same depth  #
#,silimar to get_slice in CitcomS. Note:       #
#nz_extract is within the range of [0, noz-1]  #
#to be consistent with grid setting of Paraview#
#----------------------------------------------#
def read_sphere_data(solution_cycle_init, noz, nproc, nprocz, route, nz_extract,
                     fieldtype='temperature', lnox=33, lnoy=33):
    lnoz = int((noz-1)/nprocz) + 1 ## (65-1)/2+1=33, (65-1)/4+1=17
    lnoz_extract = (nz_extract) % (lnoz-1) 
    lnno = lnox*lnoy
    nprocz_extract = int((nz_extract) / (lnoz-1)) ## 0-3, nz_extract cannot be 64, or there will be a bug !!!
    sphere_data = NP.zeros((1, 1))
    sphere_coordinate = NP.zeros((1,2)) ##needs to be from (-pi,pi) & (-pi/2,pi/2)
    for i in range(nprocz_extract, nproc, nprocz):
        vtkname ='a.proc'+str(i)+'.'+str(solution_cycle_init)+'.vts'
        vtkname =os.path.join(route,vtkname)
        field_data, coord =read_vtk_file(vtkname) ## read fields of original vts files
        extract_fd = field_data[fieldtype]
        tmp_sphere_data = extract_fd[lnoz_extract::lnoz] ##extract data from sphere: nz=nz_extract
        tmp_coordinate = coord[lnoz_extract::lnoz,:] ##extract coordinates from sphere: nz=nz_extract
        tmp_Rnorm = NP.sqrt(tmp_coordinate[:,0]**2 + tmp_coordinate[:,1]**2 + tmp_coordinate[:,2]**2)
        tmp_colati = NP.arccos(tmp_coordinate[:,2]/tmp_Rnorm)  ## (0,pi)
        tmp_longit = NP.arctan2(tmp_coordinate[:,0], tmp_coordinate[:,1])## (-pi, pi)
        tmp_Rnorm = tmp_Rnorm.reshape(lnno,1)
        tmp_colati = tmp_colati.reshape(lnno,1)
        tmp_longit = tmp_longit.reshape(lnno,1)
        tmp_sphere_data = tmp_sphere_data.reshape(lnno,1)
        tmp_sphere_coordinate = NP.hstack((tmp_longit, tmp_colati))
        sphere_data = NP.vstack((sphere_data, tmp_sphere_data)) ## the head line is useless
        sphere_coordinate = NP.vstack((sphere_coordinate, tmp_sphere_coordinate)) #the head line is useless
    return sphere_data, sphere_coordinate

def read_sphere_data_rotated(solution_cycle_init, noz, nproc, nprocz, route, nz_extract, rot_matrix, 
                             fieldtype='temperature', lnox=33, lnoy=33):
    lnoz = int((noz-1)/nprocz) + 1 ## (65-1)/2+1=33, (65-1)/4+1=17
    lnoz_extract = (nz_extract) % (lnoz-1)
    lnno = lnox*lnoy
    nprocz_extract = int((nz_extract) / (lnoz-1)) ## 0-3, nz_extract cannot be 64, or there will be a bug !!!
    sphere_data = NP.zeros((1, 1))
    sphere_coordinate = NP.zeros((1,2)) ##needs to be from (-pi,pi) & (-pi/2,pi/2)
    for i in range(nprocz_extract, nproc, nprocz):
        vtkname ='a.proc'+str(i)+'.'+str(solution_cycle_init)+'.vts'
        vtkname =os.path.join(route,vtkname)
        field_data, coord =read_vtk_file(vtkname) # read fields of original vts files
        extract_fd = field_data[fieldtype]
        tmp_sphere_data = extract_fd[lnoz_extract::lnoz]
        tmp_coordinate = coord[lnoz_extract::lnoz,:] 
        tmp_coordinate_rot = NP.matmul(tmp_coordinate, rot_matrix.T)
        tmp_Rnorm = NP.sqrt(tmp_coordinate_rot[:,0]**2 + tmp_coordinate_rot[:,1]**2 + tmp_coordinate_rot[:,2]**2)
        tmp_colati = NP.arccos(tmp_coordinate_rot[:,2]/tmp_Rnorm) ## (0,pi)
        tmp_longit = NP.arctan2(tmp_coordinate_rot[:,0], tmp_coordinate_rot[:,1])## (-pi, pi)
        tmp_Rnorm = tmp_Rnorm.reshape(lnno,1)
        tmp_colati = tmp_colati.reshape(lnno,1)
        tmp_longit = tmp_longit.reshape(lnno,1)
        tmp_sphere_data = tmp_sphere_data.reshape(lnno,1)
        tmp_sphere_coordinate = NP.hstack((tmp_longit, tmp_colati))
        sphere_data = NP.vstack((sphere_data, tmp_sphere_data)) ## the head line is useless
        sphere_coordinate = NP.vstack((sphere_coordinate, tmp_sphere_coordinate)) #the head line is useless
    return sphere_data, sphere_coordinate


#------------------------------------#
#get the velocity field at same depth#
#, silimar to get_slice in CitcomS#
#------------------------------------#
def read_sphere_velocity(solution_cycle_init, noz, nproc, nprocz, route, nz_extract,
                     fieldtype='velocity', lnox=33, lnoy=33):
    lnoz = int((noz-1)/nprocz)+1 ## (65-1)/2+1=33, (65-1)/4+1=17
    lnoz_extract = (nz_extract) % (lnoz-1)
    lnno = lnox*lnoy
    nprocz_extract = int((nz_extract) / (lnoz-1)) ## 0-3, nz_extract cannot be 64, or there will be a bug !!!
    sphere_velocity = NP.zeros((1, 3))
    sphere_coordinate = NP.zeros((1,2)) ##needs to be from (-pi,pi) & (-pi/2,pi/2)
    for i in range(nprocz_extract, nproc, nprocz):
        vtkname ='a.proc'+str(i)+'.'+str(solution_cycle_init)+'.vts'
        vtkname =os.path.join(route,vtkname)
        field_data, coord =read_vtk_file(vtkname) ## read fields of original vts files
        extract_velo = field_data[fieldtype]
        tmp_sphere_velo = extract_velo[lnoz_extract::lnoz,:] ##extract data from sphere: nz=nz_extract
        tmp_coordinate = coord[lnoz_extract::lnoz,:] ##extract coordinates from sphere: nz=nz_extract
        tmp_Rnorm = NP.sqrt(tmp_coordinate[:,0]**2 + tmp_coordinate[:,1]**2 + tmp_coordinate[:,2]**2)
        tmp_colati = NP.arccos(tmp_coordinate[:,2]/tmp_Rnorm)  ## (0,pi)
        tmp_longit = NP.arctan2(tmp_coordinate[:,0], tmp_coordinate[:,1])## (-pi, pi)
        tmp_Rnorm = tmp_Rnorm.reshape(lnno,1)
        tmp_colati = tmp_colati.reshape(lnno,1)
        tmp_longit = tmp_longit.reshape(lnno,1)
        tmp_sphere_velo = tmp_sphere_velo.reshape(lnno,3)
        tmp_sphere_coordinate = NP.hstack((tmp_longit, tmp_colati))
        sphere_velocity = NP.vstack((sphere_velocity, tmp_sphere_velo)) ## the head line is useless
        sphere_coordinate = NP.vstack((sphere_coordinate, tmp_sphere_coordinate)) #the head line is useless
    return sphere_velocity, sphere_coordinate


#------------------------------------------------------------#
#read profiles of initial condition like solidus and liquidus#
#------------------------------------------------------------#
def read_solliq_file(filename):
    with open(filename,'r') as profile:
        n=0
        radius = []
        solidus = []
        liquidus = []
        for row in profile.readlines():
            row =row.split(' ')
            radius.append(float(row[0]))
            solidus.append(float(row[1]))
            liquidus.append(float(row[2]))
            n +=1
        return radius, solidus, liquidus
    
#---read in heatflux at top---#
#Rdata usually has spare rows#
def read_topheatflux(Croute,Ccasename):
    filename = "%s/%s.qb.dat"%(Croute,Ccasename)
    fin = open(filename,'r')
    time = list()
    heatflux =list()
    sqrtvdotv =list()
    for line in fin:
        temp = line.split('  ')
        time.append(float(temp[1]))
        heatflux.append(float(temp[2]))
        sqrtvdotv.append(float(temp[3]))
    #-----drop redundant data-----#
    Rdata = NP.zeros((len(time),3))
    ftime =time[len(time)-1]+1e-9
    n_s =len(time)-1
    for i in range(len(time)-1,-1,-1):
        if time[i]<ftime:
            Rdata[n_s,0] = time[i]
            Rdata[n_s,1] = heatflux[i]
            Rdata[n_s,2] = sqrtvdotv[i]
            n_s =n_s-1
            ftime =time[i]+1e-9
    fin.close();
    return Rdata
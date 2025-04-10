#=============================================================================================
#Corrugated Cardboard to Shell V1.0 (2024/05/25) developed by ShaoKeng Liang
#This script is for corrugated cardboard homogenization in Abaqus.
#This script gets the ABDH_Matrix used in "General Shell Section" for conventional shell.
#EasyPBC Ver.1.4 developed by Sadik Lafta Omairey is used with modification in this program.
#It is recommended to edit this file on VS Code.
#Email: HenryLeung928@gmail.com
#=============================================================================================
"""
Copyright (C) [2024] [ShaoKeng Liang]

This program is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free
Software Foundation; either version 3 of the License, or any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
more details.

You should have received a copy of the GNU General Public License along
with this program; if not, see <http://www.gnu.org/licenses/>.
"""

Version='V1.0'


###
from abaqus import *
from abaqusConstants import *
from caeModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()

import numpy as np
from numpy import cos
import math
import time

start_time = time.time()


#=================parameter=====================
#User needs to modify the parameter here before running.
wallnum=2                                   #1,2,3 stands for single-wall/double-wall/triple-wall
h=[3.6,2.54]                                #height of flute, the order see <User Manual>
period=[7.7,6.41]                           #period of flute
tc=[0.221,0.221]                            #thickness of flute
tf=[0.369,0.17,0.369]                       #thickness of liner
Mat_core=(                                  #Elastic constant for liner, it is a tuple containing lists, the order is E1 E2 Nu12 G12 G13 G23
    (5740,2060,0.43,1330,133,75),
    (5710,2050,0.43,1320,133,75),
    )
Mat_face=(                                  #Elastic constant for flute, it is a tuple containing lists, the order is E1 E2 Nu12 G12 G13 G23
    (6520,2460,0.43,1550,133,75),
    (5770,2590,0.43,1500,133,75),
    (6520,2460,0.43,1550,133,75),
    )
grammage_core=[200,200]                     #grammage of flute
grammage_face=[250,250,250]                 #grammage of liner
CPUnum=6                                    #CPU numbers to use, depends on the computer configuration
cH11=1.33                                   #correction coefficient for H11 
cH22=0.94                                   #correction coefficient for H22
#==============================================


print("========== START ==========\n")
#region check before start
if wallnum>3:
    print('maximum wallnum is 3')
    error()
if len(h)!=wallnum:
    print('List size error, please check')
    error()
if len(period)!=wallnum:
    print('List size error, please check')
    error()
if len(tc)!=wallnum:
    print('List size error, please check')
    error()
if len(tf)!=wallnum+1:
    print('List size error, please check')
    error()
if len(Mat_core)!=wallnum:
    print('List size error, please check')
    error()
if len(Mat_face)!=wallnum+1:
    print('List size error, please check')
    error()
if len(grammage_core)!=wallnum:
    print('List size error, please check')
    error()
if len(grammage_face)!=wallnum+1:
    print('List size error, please check')
    error()


h = [float(x) for x in h]
period = [float(x) for x in period]
tc = [float(x) for x in tc]
tf = [float(x) for x in tf]
grammage_core = [float(x) for x in grammage_core]
grammage_face = [float(x) for x in grammage_face]
Mat_face= tuple(tuple(float(x) for x in row) for row in Mat_face)
Mat_core= tuple(tuple(float(x) for x in row) for row in Mat_core)

#endregion

Matrix={}
mass={}

#main
for j in range(wallnum):
    #region modeling 
    #part modeling
    n=1
    length=n*period[j]
    width=n*period[j]
    x=np.linspace(0, n*period[j],12*n+1)
    wave=-h[j]/2*cos(2*math.pi*x/period[j])
    peakpoint=[]
    print('===============Calculating Corrugated Flute %s===================='%(j+1))
    #
    s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
        sheetSize=200.0)
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=STANDALONE)

    for i in range(len(x)-1):
        s.Line(point1=(x[i], wave[i]), point2=(x[i+1], wave[i+1]))
    for i in range(len(x)):
        if abs(abs(wave[i])-abs(h[j]/2))<=0.01:
            peakpoint.append(i)
        else:
            continue
    for i in range(len(peakpoint)-2):
        s.Line(point1=(x[peakpoint[i]], wave[peakpoint[i]]), point2=(x[peakpoint[i+2]], wave[peakpoint[i+2]]))
    s.Line(point1=(x[0], wave[peakpoint[1]]), point2=(x[peakpoint[1]], wave[peakpoint[1]]))
    s.Line(point1=(x[-1], wave[peakpoint[-2]]), point2=(x[peakpoint[-2]], wave[peakpoint[-2]]))

    p = mdb.models['Model-1'].Part(name='Part-1', dimensionality=THREE_D, 
        type=DEFORMABLE_BODY)
    p = mdb.models['Model-1'].parts['Part-1']
    p.BaseShellExtrude(sketch=s, depth=width)
    s.unsetPrimaryObject()
    p = mdb.models['Model-1'].parts['Part-1']
    session.viewports['Viewport: 1'].setValues(displayedObject=p)
    del mdb.models['Model-1'].sketches['__profile__']

    #material property
    mdb.models['Model-1'].Material(name='Material-core2d')
    mdb.models['Model-1'].materials['Material-core2d'].Elastic(type=LAMINA, table=(
        Mat_core[j], ))
    mdb.models['Model-1'].materials['Material-core2d'].Density(table=((grammage_core[j]/tc[j]*10**-12, ), ))

    mdb.models['Model-1'].Material(name='Material-face2d')
    mdb.models['Model-1'].materials['Material-face2d'].Elastic(type=LAMINA, table=(
        Mat_face[j], ))
    mdb.models['Model-1'].materials['Material-face2d'].Density(table=((grammage_face[j]/tf[j]*10**-12, ), ))

    #section
    mdb.models['Model-1'].HomogeneousShellSection(name='Section-face', 
        preIntegrate=OFF, material='Material-face2d', thicknessType=UNIFORM, 
        thickness=tf[j], thicknessField='', nodalThicknessField='', 
        idealization=NO_IDEALIZATION, poissonDefinition=DEFAULT, 
        thicknessModulus=None, temperature=GRADIENT, useDensity=OFF, 
        integrationRule=SIMPSON, numIntPts=5)
    mdb.models['Model-1'].HomogeneousShellSection(name='Section-core', 
        preIntegrate=OFF, material='Material-core2d', thicknessType=UNIFORM, 
        thickness=tc[j], thicknessField='', nodalThicknessField='', 
        idealization=NO_IDEALIZATION, poissonDefinition=DEFAULT, 
        thicknessModulus=None, temperature=GRADIENT, useDensity=OFF, 
        integrationRule=SIMPSON, numIntPts=5)    
    #assign    
    p = mdb.models['Model-1'].parts['Part-1']
    f = p.faces
    faces = f.getSequenceFromMask(mask=('[#1c00 ]', ), )
    region = regionToolset.Region(faces=faces)
    p = mdb.models['Model-1'].parts['Part-1']
    p.SectionAssignment(region=region, sectionName='Section-face', offset=0.0, 
        offsetType=MIDDLE_SURFACE, offsetField='', 
        thicknessAssignment=FROM_SECTION)    
        
    p = mdb.models['Model-1'].parts['Part-1']
    f = p.faces
    faces = f.getSequenceFromMask(mask=('[#3ff ]', ), )
    region = regionToolset.Region(faces=faces)
    p = mdb.models['Model-1'].parts['Part-1']
    p.SectionAssignment(region=region, sectionName='Section-core', offset=0.0, 
        offsetType=MIDDLE_SURFACE, offsetField='', 
        thicknessAssignment=FROM_SECTION)    

    #assembly
    a1 = mdb.models['Model-1'].rootAssembly
    a1.DatumCsysByDefault(CARTESIAN)
    p = mdb.models['Model-1'].parts['Part-1']

    a1.Instance(name='Part-1-1', part=p, dependent=ON)
    a1 = mdb.models['Model-1'].rootAssembly
    a1.rotate(instanceList=('Part-1-1', ), axisPoint=(0.0, 0.0, 0.0), 
        axisDirection=(10.0, 0.0, 0.0), angle=90.0)

    a1 = mdb.models['Model-1'].rootAssembly
    a1.translate(instanceList=('Part-1-1', ), vector=(0.0, period[j], h[j]/2))

    #mesh
    p = mdb.models['Model-1'].parts['Part-1']
    e = p.edges
    pickedEdges = e.getSequenceFromMask(mask=('[#0 #120 ]', ), )
    p.seedEdgeByNumber(edges=pickedEdges, number=6, constraint=FINER)
    pickedEdges = e.getSequenceFromMask(mask=('[#0 #8 ]', ), )
    p.seedEdgeByNumber(edges=pickedEdges, number=12, constraint=FINER)
    p = mdb.models['Model-1'].parts['Part-1']
    p.generateMesh()

    #orientation
    p = mdb.models['Model-1'].parts['Part-1']
    f = p.faces
    faces = f.getSequenceFromMask(mask=('[#1c00 ]', ), )
    region = regionToolset.Region(faces=faces)
    orientation=None
    mdb.models['Model-1'].parts['Part-1'].MaterialOrientation(region=region, 
        orientationType=GLOBAL, axis=AXIS_3, additionalRotationType=ROTATION_NONE, 
        localCsys=None, fieldName='')

    p = mdb.models['Model-1'].parts['Part-1']
    f = p.faces
    faces = f.getSequenceFromMask(mask=('[#3ff ]', ), )
    region = regionToolset.Region(faces=faces)
    p = mdb.models['Model-1'].parts['Part-1']
    s = p.faces
    side1Faces = s.getSequenceFromMask(mask=('[#3ff ]', ), )
    normalAxisRegion = regionToolset.Region(side1Faces=side1Faces)
    p = mdb.models['Model-1'].parts['Part-1']
    e = p.edges
    edges = e.getSequenceFromMask(mask=('[#31249891 #1 ]', ), )
    primaryAxisRegion = regionToolset.Region(edges=edges)
    mdb.models['Model-1'].parts['Part-1'].MaterialOrientation(region=region, 
        orientationType=DISCRETE, axis=AXIS_3, normalAxisDefinition=SURFACE, 
        normalAxisRegion=normalAxisRegion, flipNormalDirection=False, 
        normalAxisDirection=AXIS_3, primaryAxisDefinition=EDGE, 
        primaryAxisRegion=primaryAxisRegion, primaryAxisDirection=AXIS_1, 
        flipPrimaryDirection=False, additionalRotationType=ROTATION_NONE, 
        angle=0.0, additionalRotationField='')
        
    #finding specific nodes
    a = mdb.models['Model-1'].rootAssembly
    a.regenerate()
    Nodeset = mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].nodes
    if len(Nodeset)==0:
        error()
    dict_nodes ={}
    tor=1E-3
    for i in Nodeset:
        if abs(i.coordinates[0]-0)<tor:
            if abs(i.coordinates[1]-0)<tor:
                if abs(i.coordinates[2]-0)<tor:
                        dict_nodes[7] =i.label
                elif abs(i.coordinates[2]-h[j])<tor:
                        dict_nodes[5] =i.label
                else:
                    continue
            elif abs(i.coordinates[1]-width/2)<tor:
                if abs(i.coordinates[2]-0)<tor:
                    dict_nodes[8] =i.label
                elif abs(i.coordinates[2]-h[j])<tor:
                    dict_nodes[6] =i.label
                else:
                    continue
            else:
                continue
        elif abs(i.coordinates[0]-length/2)<tor:
            if abs(i.coordinates[1]-0)<tor:
                if abs(i.coordinates[2]-0)<tor:
                        dict_nodes[3] =i.label
                elif abs(i.coordinates[2]-h[j])<tor:
                        dict_nodes[1] =i.label
                else:
                    continue
            elif abs(i.coordinates[1]-width/2)<tor:
                if abs(i.coordinates[2]-0)<tor:
                    dict_nodes[4] =i.label
                elif abs(i.coordinates[2]-h[j])<tor:
                    dict_nodes[2] =i.label
                else:
                    continue   
            else:
                continue 
        else:
            continue     
    
    sorted_keys = sorted(dict_nodes.keys())
    interested_nodes = [dict_nodes[key] for key in sorted_keys]



    #: Model-2 for bending test
    mdb.Model(name='Model-2', objectToCopy=mdb.models['Model-1'])

    #endregion
 
    #region EasyPBC
    
    #region EasyPBC definition      
    def feasypbc(part,inst,meshsens,E33,CPU):
        import os
        path = os.getcwd()
        start = time.time()

        modelName = part
        instanceName = inst
        upperName= inst.upper()
        CPUs = int(round(CPU))
        
        a = mdb.models[modelName].rootAssembly      
        Nodeset = mdb.models[modelName].rootAssembly.instances[instanceName].nodes

        ## Start of sets creation ## 
        #region variable               
        j = 0
        x=[]
        y=[]
        z=[]
        c1=[]
        c2=[]
        c3=[]
        c4=[]
        c5=[]
        c6=[]
        c7=[]
        c8=[]
        Max=[]
        ftedgexyz={}
        btedgexyz={}
        fbedgexyz={}
        bbedgexyz={}
        fledgexyz={}
        bledgexyz={}
        fredgexyz={}
        bredgexyz={}
        ltedgexyz={}
        rtedgexyz={}
        lbedgexyz={}
        rbedgexyz={}
        frontsxyz={}
        backsxyz={}
        topsxyz={}
        botsxyz={}
        leftsxyz={}
        rightsxyz={}
        frontbcxyz={}
        backbcxyz={}
        topbcxyz={}
        botbcxyz={}
        leftbcxyz={}
        rightbcxyz={}
        ftedge=[]
        btedge=[]
        fbedge=[]
        bbedge=[]
        fledge=[]
        fredge=[]
        bledge=[]
        bredge=[]
        ltedge=[]
        lbedge=[]
        rtedge=[]
        rbedge=[]
        fronts=[]
        backs=[]
        lefts=[]
        rights=[]
        tops=[]
        bots=[]
        backs=[]
        frontbc=[]
        backbc=[]
        leftbc=[]
        rightbc=[]
        topbc=[]
        botbc=[]
        backbc=[]
        errorset=[]
        coc1={}
        coc2={}
        coc3={}
        coc4={}
        coc5={}
        coc6={}
        coc7={}
        coc8={}

        #endregion

        print ('-------- Start of Modified_EasyPBC --------')

        ## Identifying RVE size ##    
        for i in Nodeset:
            x.insert(j,i.coordinates[0])
            y.insert(j,i.coordinates[1])
            z.insert(j,i.coordinates[2])
            j=j+1

        Max = max(x)
        May = max(y)
        Maz = max(z)
        Mnx = min(x)
        Mny = min(y)
        Mnz = min(z)

        ## 3D model ##
        
        L=abs(Max-Mnx)
        H=abs(May-Mny)
        W=abs(Maz-Mnz)
        
        Dispx = L*0.2
        Dispy = H*0.2
        Dispz = W*0.2
        
        ## Creating Ref. Points ##
        for i in a.features.keys():
            if i.startswith('RP'):
                del a.features['%s' % (i)]
        a.ReferencePoint(point=(Max+0.8*abs(Max-Mnx), May-0.5*(May-Mny), Maz-0.5*(Maz-Mnz)))  ## RP6: G23
        a.ReferencePoint(point=(Max+0.6*abs(Max-Mnx), May-0.5*(May-Mny), Maz-0.5*(Maz-Mnz)))  ## RP5: G13
        a.ReferencePoint(point=(Max+0.4*abs(Max-Mnx), May-0.5*(May-Mny), Maz-0.5*(Maz-Mnz)))  ## RP4: G12
        a.ReferencePoint(point=(Max+0.2*abs(Max-Mnx), May-0.5*(May-Mny), Maz-0.5*(Maz-Mnz)))  ## RP3: Rigid body movement X-axis
        a.ReferencePoint(point=(Max-0.5*(Max-Mnx), May-0.5*(May-Mny), Maz+0.2*abs(Maz-Mnz)))  ## RP2: Rigid body movement Z-axis
        a.ReferencePoint(point=(Max-0.5*(Max-Mnx), May+0.2*abs(May-Mny), Maz-0.5*(Maz-Mnz)))  ## RP1: Rigid body movement Y-axis

        r1 = a.referencePoints

        ## Naming Ref. Points ##
        d=1
        for i in r1.keys():
            refPoints1=(r1[i], )
            a.Set(referencePoints=refPoints1, name='RP%s' % (d))
            d=d+1
            
        ## Identifying boundary nodes ##
        for i in Nodeset:
            if (Mnx+meshsens) < i.coordinates[0] < (Max-meshsens) and (Mny+meshsens) < i.coordinates[1] < (May-meshsens) and (Mnz+meshsens) < i.coordinates[2] < (Maz-meshsens):
                continue
            if abs(i.coordinates[0]-Max)<=meshsens:
                frontbcxyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]	
            if abs(i.coordinates[0]-Mnx)<=meshsens:
                backbcxyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]	
            if abs(i.coordinates[2]-Maz)<=meshsens:
                leftbcxyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]
            if abs(i.coordinates[2]-Mnz)<=meshsens:
                rightbcxyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]	
            if abs(i.coordinates[1]-May)<=meshsens:
                topbcxyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]	
            if abs(i.coordinates[1]-Mny)<=meshsens:
                botbcxyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]
            if abs(i.coordinates[0]-Max)<=meshsens and abs(i.coordinates[1]-May)<=meshsens and abs(i.coordinates[2]-Maz)<=meshsens:
                c1.insert(0,i.label)
                coc1[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]
            if abs(i.coordinates[0]-Mnx)<=meshsens and abs(i.coordinates[1]-May)<=meshsens and abs(i.coordinates[2]-Maz)<=meshsens:
                c2.insert(0,i.label)
                coc2[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]
            if abs(i.coordinates[0]-Mnx)<=meshsens and abs(i.coordinates[1]-May)<=meshsens and abs(i.coordinates[2]-Mnz)<=meshsens:
                c3.insert(0,i.label)
                coc3[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]
            if abs(i.coordinates[0]-Max)<=meshsens and abs(i.coordinates[1]-May)<=meshsens and abs(i.coordinates[2]-Mnz)<=meshsens:
                c4.insert(0,i.label)
                coc4[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]
            if abs(i.coordinates[0]-Max)<=meshsens and abs(i.coordinates[1]-Mny)<=meshsens and abs(i.coordinates[2]-Maz)<=meshsens:
                c5.insert(0,i.label)
                coc5[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]
            if abs(i.coordinates[0]-Mnx)<=meshsens and abs(i.coordinates[1]-Mny)<=meshsens and abs(i.coordinates[2]-Maz)<=meshsens:
                c6.insert(0,i.label)
                coc6[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]
            if abs(i.coordinates[0]-Mnx)<=meshsens and abs(i.coordinates[1]-Mny)<=meshsens and abs(i.coordinates[2]-Mnz)<=meshsens:
                c7.insert(0,i.label)
                coc7[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]
            if abs(i.coordinates[0]-Max)<=meshsens and abs(i.coordinates[1]-Mny)<=meshsens and abs(i.coordinates[2]-Mnz)<=meshsens:
                c8.insert(0,i.label)
                coc8[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]
            if abs(i.coordinates[0]-Max)<=meshsens and abs(i.coordinates[1]-May)<=meshsens and abs(i.coordinates[2]-Maz)>meshsens and abs(i.coordinates[2]-Mnz)>meshsens:
                ftedgexyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]
            if abs(i.coordinates[0]-Max)<=meshsens and abs(i.coordinates[1]-Mny)<=meshsens and abs(i.coordinates[2]-Maz)>meshsens and abs(i.coordinates[2]-Mnz)>meshsens:
                fbedgexyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]
            if abs(i.coordinates[0]-Mnx)<=meshsens and abs(i.coordinates[1]-May)<=meshsens and abs(i.coordinates[2]-Maz)>meshsens and abs(i.coordinates[2]-Mnz)>meshsens:
                btedgexyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]
            if abs(i.coordinates[0]-Mnx)<=meshsens and abs(i.coordinates[1]-Mny)<=meshsens and abs(i.coordinates[2]-Maz)>meshsens and abs(i.coordinates[2]-Mnz)>meshsens:
                bbedgexyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]	
            if abs(i.coordinates[0]-Max)<=meshsens and abs(i.coordinates[2]-Maz)<=meshsens and abs(i.coordinates[1]-May)>meshsens and abs(i.coordinates[1]-Mny)>meshsens:
                fledgexyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]	
            if abs(i.coordinates[0]-Max)<=meshsens and abs(i.coordinates[2]-Mnz)<=meshsens and abs(i.coordinates[1]-May)>meshsens and abs(i.coordinates[1]-Mny)>meshsens:
                fredgexyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]	
            if abs(i.coordinates[0]-Mnx)<=meshsens and abs(i.coordinates[2]-Maz)<=meshsens and abs(i.coordinates[1]-May)>meshsens and abs(i.coordinates[1]-Mny)>meshsens:
                bledgexyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]	
            if abs(i.coordinates[0]-Mnx)<=meshsens and abs(i.coordinates[2]-Mnz)<=meshsens and abs(i.coordinates[1]-May)>meshsens and abs(i.coordinates[1]-Mny)>meshsens:
                bredgexyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]	
            if abs(i.coordinates[2]-Maz)<=meshsens and abs(i.coordinates[1]-May)<=meshsens and abs(i.coordinates[0]-Max)>meshsens and abs(i.coordinates[0]-Mnx)>meshsens:
                ltedgexyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]	
            if abs(i.coordinates[2]-Maz)<=meshsens and abs(i.coordinates[1]-Mny)<=meshsens and abs(i.coordinates[0]-Max)>meshsens and abs(i.coordinates[0]-Mnx)>meshsens:
                lbedgexyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]	
            if abs(i.coordinates[2]-Mnz)<=meshsens and abs(i.coordinates[1]-May)<=meshsens and abs(i.coordinates[0]-Max)>meshsens and abs(i.coordinates[0]-Mnx)>meshsens:
                rtedgexyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]	
            if abs(i.coordinates[2]-Mnz)<=meshsens and abs(i.coordinates[1]-Mny)<=meshsens and abs(i.coordinates[0]-Max)>meshsens and abs(i.coordinates[0]-Mnx)>meshsens:
                rbedgexyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]	
            if abs(i.coordinates[0]-Max)<=meshsens and abs(i.coordinates[1]-May)>meshsens and abs(i.coordinates[1]-Mny)>meshsens and abs(i.coordinates[2]-Maz)>meshsens and abs(i.coordinates[2]-Mnz)>meshsens:
                frontsxyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]
            if abs(i.coordinates[0]-Mnx)<=meshsens and abs(i.coordinates[1]-May)>meshsens and abs(i.coordinates[1]-Mny)>meshsens and abs(i.coordinates[2]-Maz)>meshsens and abs(i.coordinates[2]-Mnz)>meshsens:
                backsxyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]] 
            if abs(i.coordinates[2]-Maz)<=meshsens and abs(i.coordinates[1]-May)>meshsens and abs(i.coordinates[1]-Mny)>meshsens and abs(i.coordinates[0]-Max)>meshsens and abs(i.coordinates[0]-Mnx)>meshsens:
                leftsxyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]	
            if abs(i.coordinates[2]-Mnz)<=meshsens and abs(i.coordinates[1]-May)>meshsens and abs(i.coordinates[1]-Mny)>meshsens and abs(i.coordinates[0]-Max)>meshsens and abs(i.coordinates[0]-Mnx)>meshsens:
                rightsxyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]   
            if abs(i.coordinates[1]-May)<=meshsens and abs(i.coordinates[0]-Max)>meshsens and abs(i.coordinates[0]-Mnx)>meshsens and abs(i.coordinates[2]-Maz)>meshsens and abs(i.coordinates[2]-Mnz)>meshsens:
                topsxyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]
            if abs(i.coordinates[1]-Mny)<=meshsens and abs(i.coordinates[0]-Max)>meshsens and abs(i.coordinates[0]-Mnx)>meshsens and abs(i.coordinates[2]-Maz)>meshsens and abs(i.coordinates[2]-Mnz)>meshsens:
                botsxyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]

        #region Sorting and appending sets ##
        for i in frontsxyz.keys():
                for k in backsxyz.keys():
                        if abs(frontsxyz[i][1] - backsxyz[k][1])<=meshsens and abs(frontsxyz[i][2] - backsxyz[k][2])<=meshsens:
                                fronts.append(i)
                                backs.append(k)

        for i in topsxyz.keys():
            for k in botsxyz.keys():
                if abs(topsxyz[i][0] - botsxyz[k][0]) <=meshsens and abs(topsxyz[i][2] - botsxyz[k][2]) <=meshsens:
                    tops.append(i)
                    bots.append(k)

        for i in leftsxyz.keys():
            for k in rightsxyz.keys():
                if abs(leftsxyz[i][0] - rightsxyz[k][0])<=meshsens and abs(leftsxyz[i][1] - rightsxyz[k][1]) <=meshsens:
                    lefts.append(i)
                    rights.append(k)

        for i in frontbcxyz.keys():
            for k in backbcxyz.keys():
                if abs(frontbcxyz[i][1] - backbcxyz[k][1])<=meshsens and abs(frontbcxyz[i][2] - backbcxyz[k][2])<=meshsens:
                    frontbc.append(i)
                    backbc.append(k)

        for i in topbcxyz.keys():
            for k in botbcxyz.keys():
                if abs(topbcxyz[i][0] - botbcxyz[k][0]) <=meshsens and abs(topbcxyz[i][2] - botbcxyz[k][2]) <=meshsens:
                    topbc.append(i)
                    botbc.append(k)

        for i in leftbcxyz.keys():
            for k in rightbcxyz.keys():
                if abs(leftbcxyz[i][0] - rightbcxyz[k][0])<=meshsens and abs(leftbcxyz[i][1] - rightbcxyz[k][1]) <=meshsens:
                    leftbc.append(i)
                    rightbc.append(k)

        for i in ftedgexyz.keys():
            for k in btedgexyz.keys():
                if abs(ftedgexyz[i][1] - btedgexyz[k][1])<=meshsens and abs(ftedgexyz[i][2] - btedgexyz[k][2])<=meshsens:
                    ftedge.append(i)
                    btedge.append(k)
        for i in btedge:
            for k in bbedgexyz.keys():
                if abs(btedgexyz[i][0] - bbedgexyz[k][0]) <=meshsens and abs(btedgexyz[i][2] - bbedgexyz[k][2]) <=meshsens:
                    bbedge.append(k)    
        for i in bbedge:
            for k in fbedgexyz.keys():
                if abs(bbedgexyz[i][1] - fbedgexyz[k][1]) <=meshsens and abs(bbedgexyz[i][2] - fbedgexyz[k][2]) <=meshsens:
                    fbedge.append(k) 

        for i in ltedgexyz.keys():
            for k in rtedgexyz.keys():
                if abs(ltedgexyz[i][0] - rtedgexyz[k][0])<=meshsens and abs(ltedgexyz[i][1] - rtedgexyz[k][1])<=meshsens:
                    ltedge.append(i)
                    rtedge.append(k)
        for i in rtedge:
            for k in rbedgexyz.keys():
                if abs(rtedgexyz[i][0] - rbedgexyz[k][0])<=meshsens and abs(rtedgexyz[i][2] - rbedgexyz[k][2])<=meshsens:
                    rbedge.append(k)    
        for i in rbedge:
            for k in lbedgexyz.keys():
                if abs(rbedgexyz[i][0] - lbedgexyz[k][0])<=meshsens and abs(rbedgexyz[i][1] - lbedgexyz[k][1])<=meshsens:
                    lbedge.append(k) 

        for i in fledgexyz.keys():
            for k in bledgexyz.keys():
                if abs(fledgexyz[i][1] - bledgexyz[k][1])<=meshsens and abs(fledgexyz[i][2] - bledgexyz[k][2])<=meshsens:
                    fledge.append(i)
                    bledge.append(k)
        for i in bledge:
            for k in bredgexyz.keys():
                if abs(bledgexyz[i][0] - bredgexyz[k][0])<=meshsens and abs(bledgexyz[i][1] - bredgexyz[k][1])<=meshsens:
                    bredge.append(k)    
        for i in bredge:
            for k in fredgexyz.keys():
                if abs(bredgexyz[i][1] - fredgexyz[k][1])<=meshsens and abs(bredgexyz[i][2] - fredgexyz[k][2])<=meshsens:
                    fredge.append(k) 
        #endregion

        #region Creating ABAQUS sets ##
        a.SetFromNodeLabels(name='c1', nodeLabels=((instanceName,c1),))
        a.SetFromNodeLabels(name='c2', nodeLabels=((instanceName,c2),))
        a.SetFromNodeLabels(name='c3', nodeLabels=((instanceName,c3),))
        a.SetFromNodeLabels(name='c4', nodeLabels=((instanceName,c4),))
        a.SetFromNodeLabels(name='c5', nodeLabels=((instanceName,c5),))
        a.SetFromNodeLabels(name='c6', nodeLabels=((instanceName,c6),))
        a.SetFromNodeLabels(name='c7', nodeLabels=((instanceName,c7),))
        a.SetFromNodeLabels(name='c8', nodeLabels=((instanceName,c8),))
        #a.SetFromNodeLabels(name='ftedge', nodeLabels=((instanceName,ftedge),))#These codes are commented out because corrugated cardboard model does not have such nodes.
        #a.SetFromNodeLabels(name='fbedge', nodeLabels=((instanceName,fbedge),))
        #a.SetFromNodeLabels(name='btedge', nodeLabels=((instanceName,btedge),))
        #a.SetFromNodeLabels(name='bbedge', nodeLabels=((instanceName,bbedge),))
        a.SetFromNodeLabels(name='fledge', nodeLabels=((instanceName,fledge),))
        a.SetFromNodeLabels(name='fredge', nodeLabels=((instanceName,fredge),))
        a.SetFromNodeLabels(name='bledge', nodeLabels=((instanceName,bledge),))
        a.SetFromNodeLabels(name='bredge', nodeLabels=((instanceName,bredge),))
        a.SetFromNodeLabels(name='ltedge', nodeLabels=((instanceName,ltedge),))
        a.SetFromNodeLabels(name='lbedge', nodeLabels=((instanceName,lbedge),))
        a.SetFromNodeLabels(name='rtedge', nodeLabels=((instanceName,rtedge),))
        a.SetFromNodeLabels(name='rbedge', nodeLabels=((instanceName,rbedge),))
        #a.SetFromNodeLabels(name='fronts', nodeLabels=((instanceName,fronts),))
        #a.SetFromNodeLabels(name='backs', nodeLabels=((instanceName,backs),))
        a.SetFromNodeLabels(name='lefts', nodeLabels=((instanceName,lefts),))
        a.SetFromNodeLabels(name='rights', nodeLabels=((instanceName,rights),))
        a.SetFromNodeLabels(name='tops', nodeLabels=((instanceName,tops),))
        a.SetFromNodeLabels(name='bots', nodeLabels=((instanceName,bots),))
        a.SetFromNodeLabels(name='frontbc', nodeLabels=((instanceName,frontbc),))
        a.SetFromNodeLabels(name='backbc', nodeLabels=((instanceName,backbc),))
        a.SetFromNodeLabels(name='leftbc', nodeLabels=((instanceName,leftbc),))
        a.SetFromNodeLabels(name='rightbc', nodeLabels=((instanceName,rightbc),))
        a.SetFromNodeLabels(name='topbc', nodeLabels=((instanceName,topbc),))
        a.SetFromNodeLabels(name='botbc', nodeLabels=((instanceName,botbc),))
        
        #endregion

        #Step
        mdb.models[modelName].StaticStep(name='Step-1', previous='Initial')

        #region Creating single-node ABAQUS sets ##                
        for i,k in zip(tops,bots):
            a.SetFromNodeLabels(name='tops%s' % (i), nodeLabels=((instanceName,[i]),))
            a.SetFromNodeLabels(name='bots%s' % (k), nodeLabels=((instanceName,[k]),))

        for i,k in zip(fronts,backs):
            a.SetFromNodeLabels(name='fronts%s' % (i), nodeLabels=((instanceName,[i]),))
            a.SetFromNodeLabels(name='backs%s' % (k), nodeLabels=((instanceName,[k]),))

        for i,k in zip(lefts,rights):
            a.SetFromNodeLabels(name='lefts%s' % (i), nodeLabels=((instanceName,[i]),))
            a.SetFromNodeLabels(name='rights%s' % (k), nodeLabels=((instanceName,[k]),))

        for i,k,j,l in zip(ftedge,btedge,bbedge,fbedge):
            a.SetFromNodeLabels(name='ftedge%s' % (i), nodeLabels=((instanceName,[i]),))
            a.SetFromNodeLabels(name='btedge%s' % (k), nodeLabels=((instanceName,[k]),))
            a.SetFromNodeLabels(name='bbedge%s' % (j), nodeLabels=((instanceName,[j]),))
            a.SetFromNodeLabels(name='fbedge%s' % (l), nodeLabels=((instanceName,[l]),))

        for i,k,j,l in zip(fledge,bledge,bredge,fredge):
            a.SetFromNodeLabels(name='fledge%s' % (i), nodeLabels=((instanceName,[i]),))
            a.SetFromNodeLabels(name='bledge%s' % (k), nodeLabels=((instanceName,[k]),))
            a.SetFromNodeLabels(name='bredge%s' % (j), nodeLabels=((instanceName,[j]),))
            a.SetFromNodeLabels(name='fredge%s' % (l), nodeLabels=((instanceName,[l]),))

        for i,k,j,l in zip(ltedge,lbedge,rbedge,rtedge):
            a.SetFromNodeLabels(name='ltedge%s' % (i), nodeLabels=((instanceName,[i]),))
            a.SetFromNodeLabels(name='lbedge%s' % (k), nodeLabels=((instanceName,[k]),))
            a.SetFromNodeLabels(name='rbedge%s' % (j), nodeLabels=((instanceName,[j]),))
            a.SetFromNodeLabels(name='rtedge%s' % (l), nodeLabels=((instanceName,[l]),))

        for i,k in zip(topbc,botbc):
            a.SetFromNodeLabels(name='topbc%s' % (i), nodeLabels=((instanceName,[i]),))
            a.SetFromNodeLabels(name='botbc%s' % (k), nodeLabels=((instanceName,[k]),))

        for i,k in zip(frontbc,backbc):
            a.SetFromNodeLabels(name='frontbc%s' % (i), nodeLabels=((instanceName,[i]),))
            a.SetFromNodeLabels(name='backbc%s' % (k), nodeLabels=((instanceName,[k]),))

        for i,k in zip(leftbc,rightbc):
            a.SetFromNodeLabels(name='leftbc%s' % (i), nodeLabels=((instanceName,[i]),))
            a.SetFromNodeLabels(name='rightbc%s' % (k), nodeLabels=((instanceName,[k]),))
        #endregion

        #region Creating constraints for elastic moduli ##
            
        for i in mdb.models[modelName].constraints.keys():
                del mdb.models[modelName].constraints[i]
                
        for i,k in zip(tops,bots):
            mdb.models[modelName].Equation(name='E-1-tops-bots%s'%i, terms=((1.0, 'tops%s'%i, 1), (-1.0, 'bots%s'%k, 1)))
        for i,k in zip(tops,bots):
            mdb.models[modelName].Equation(name='E-2-tops-bots%s'%i, terms=((1.0, 'tops%s'%i, 2), (-1.0, 'bots%s'%k, 2),(-1.0, 'RP5', 2)))
        for i,k in zip(tops,bots):
            mdb.models[modelName].Equation(name='E-3-tops-bots%s'%i, terms=((1.0, 'tops%s'%i, 3), (-1.0, 'bots%s'%k, 3)))

        for i,k in zip(lefts,rights):
            mdb.models[modelName].Equation(name='E-1-lefts-rights%s'%i, terms=((1.0, 'lefts%s'%i, 1), (-1.0, 'rights%s'%k, 1)))
        for i,k in zip(lefts,rights):
            mdb.models[modelName].Equation(name='E-2-lefts-rights%s'%i, terms=((1.0, 'lefts%s'%i, 2), (-1.0, 'rights%s'%k, 2)))
        for i,k in zip(lefts,rights):
            mdb.models[modelName].Equation(name='E-3-lefts-rights%s'%i, terms=((1.0, 'lefts%s'%i, 3), (-1.0, 'rights%s'%k, 3),(-1.0, 'RP6', 3)))

        for i,k in zip(fronts,backs):
            mdb.models[modelName].Equation(name='E-1-fronts-backs%s'%i, terms=((1.0, 'fronts%s'%i, 1), (-1.0, 'backs%s'%k, 1),(-1.0, 'RP4', 1)))
        for i,k in zip(fronts,backs):
            mdb.models[modelName].Equation(name='E-2-fronts-backs%s'%i, terms=((1.0, 'fronts%s'%i, 2), (-1.0, 'backs%s'%k, 2)))
        for i,k in zip(fronts,backs):
            mdb.models[modelName].Equation(name='E-3-fronts-backs%s'%i, terms=((1.0, 'fronts%s'%i, 3), (-1.0, 'backs%s'%k, 3)))


        mdb.models[modelName].Equation(name='E-1-c12', terms=((1.0, 'c6', 1), (-1.0, 'c2', 1)))
        mdb.models[modelName].Equation(name='E-1-c23', terms=((1.0, 'c2', 1), (-1.0, 'c3', 1)))
        mdb.models[modelName].Equation(name='E-1-c34', terms=((1.0, 'c3', 1), (-1.0, 'c4', 1),(1.0, 'RP4', 1)))
        mdb.models[modelName].Equation(name='E-1-c45', terms=((1.0, 'c4', 1), (-1.0, 'c8', 1)))
        mdb.models[modelName].Equation(name='E-1-c56', terms=((1.0, 'c8', 1), (-1.0, 'c5', 1)))
        mdb.models[modelName].Equation(name='E-1-c67', terms=((1.0, 'c5', 1), (-1.0, 'c1', 1)))
        mdb.models[modelName].Equation(name='E-1-c78', terms=((1.0, 'c1', 1), (-1.0, 'c7', 1),(-1.0, 'RP4', 1)))
        
        mdb.models[modelName].Equation(name='E-2-c12', terms=((1.0, 'c6', 2), (-1.0, 'c2', 2),(1.0, 'RP5', 2)))
        mdb.models[modelName].Equation(name='E-2-c23', terms=((1.0, 'c2', 2), (-1.0, 'c3', 2)))
        mdb.models[modelName].Equation(name='E-2-c34', terms=((1.0, 'c3', 2), (-1.0, 'c4', 2)))
        mdb.models[modelName].Equation(name='E-2-c45', terms=((1.0, 'c4', 2), (-1.0, 'c8', 2),(-1.0, 'RP5', 2)))
        mdb.models[modelName].Equation(name='E-2-c56', terms=((1.0, 'c8', 2), (-1.0, 'c5', 2)))
        mdb.models[modelName].Equation(name='E-2-c67', terms=((1.0, 'c5', 2), (-1.0, 'c1', 2),(1.0, 'RP5', 2)))
        mdb.models[modelName].Equation(name='E-2-c78', terms=((1.0, 'c1', 2), (-1.0, 'c7', 2),(-1.0, 'RP5', 2)))

        mdb.models[modelName].Equation(name='E-3-c12', terms=((1.0, 'c6', 3), (-1.0, 'c2', 3)))
        mdb.models[modelName].Equation(name='E-3-c23', terms=((1.0, 'c2', 3), (-1.0, 'c3', 3),(-1.0, 'RP6', 3)))
        mdb.models[modelName].Equation(name='E-3-c34', terms=((1.0, 'c3', 3), (-1.0, 'c4', 3)))
        mdb.models[modelName].Equation(name='E-3-c45', terms=((1.0, 'c4', 3), (-1.0, 'c8', 3)))
        mdb.models[modelName].Equation(name='E-3-c56', terms=((1.0, 'c8', 3), (-1.0, 'c5', 3),(1.0, 'RP6', 3)))
        mdb.models[modelName].Equation(name='E-3-c67', terms=((1.0, 'c5', 3), (-1.0, 'c1', 3)))
        mdb.models[modelName].Equation(name='E-3-c78', terms=((1.0, 'c1', 3), (-1.0, 'c7', 3),(-1.0, 'RP6', 3)))
                

        for i,k,j,l in zip(ftedge,btedge,bbedge,fbedge):
            mdb.models[modelName].Equation(name='E-1-ftedge-btedge%s'%i, terms=((1.0, 'ftedge%s'%i, 1), (-1.0, 'btedge%s'%k, 1),(-1.0, 'RP4', 1)))
            mdb.models[modelName].Equation(name='E-1-btedge-bbedge%s'%k, terms=((1.0, 'btedge%s'%k, 1), (-1.0, 'bbedge%s'%j, 1)))
            mdb.models[modelName].Equation(name='E-1-bbedge-fbedge%s'%j, terms=((1.0, 'bbedge%s'%j, 1), (-1.0, 'fbedge%s'%l, 1),(1.0, 'RP4', 1)))
        for i,k,j,l in zip(ftedge,btedge,bbedge,fbedge):
            mdb.models[modelName].Equation(name='E-2-ftedge-btedge%s'%i, terms=((1.0, 'ftedge%s'%i, 2), (-1.0, 'btedge%s'%k, 2)))
            mdb.models[modelName].Equation(name='E-2-btedge-bbedge%s'%k, terms=((1.0, 'btedge%s'%k, 2), (-1.0, 'bbedge%s'%j, 2),(-1.0, 'RP5', 2)))
            mdb.models[modelName].Equation(name='E-2-bbedge-fbedge%s'%j, terms=((1.0, 'bbedge%s'%j, 2), (-1.0, 'fbedge%s'%l, 2)))
        for i,k,j,l in zip(ftedge,btedge,bbedge,fbedge):
            mdb.models[modelName].Equation(name='E-3-ftedge-btedge%s'%i, terms=((1.0, 'ftedge%s'%i, 3), (-1.0, 'btedge%s'%k, 3)))
            mdb.models[modelName].Equation(name='E-3-btedge-bbedge%s'%k, terms=((1.0, 'btedge%s'%k, 3), (-1.0, 'bbedge%s'%j, 3)))
            mdb.models[modelName].Equation(name='E-3-bbedge-fbedge%s'%j, terms=((1.0, 'bbedge%s'%j, 3), (-1.0, 'fbedge%s'%l, 3)))

        for i,k,j,l in zip(fledge,bledge,bredge,fredge):
            mdb.models[modelName].Equation(name='E-1-fledge-bledge%s'%i, terms=((1.0, 'fledge%s'%i, 1), (-1.0, 'bledge%s'%k, 1),(-1.0, 'RP4', 1)))
            mdb.models[modelName].Equation(name='E-1-bledge-bredge%s'%k, terms=((1.0, 'bledge%s'%k, 1), (-1.0, 'bredge%s'%j, 1)))
            mdb.models[modelName].Equation(name='E-1-bredge-fredge%s'%j, terms=((1.0, 'bredge%s'%j, 1), (-1.0, 'fredge%s'%l, 1),(1.0, 'RP4', 1)))
        for i,k,j,l in zip(fledge,bledge,bredge,fredge):
            mdb.models[modelName].Equation(name='E-2-fledge-bledge%s'%i, terms=((1.0, 'fledge%s'%i, 2), (-1.0, 'bledge%s'%k, 2)))
            mdb.models[modelName].Equation(name='E-2-bledge-bredge%s'%k, terms=((1.0, 'bledge%s'%k, 2), (-1.0, 'bredge%s'%j, 2)))
            mdb.models[modelName].Equation(name='E-2-bredge-fredge%s'%j, terms=((1.0, 'bredge%s'%j, 2), (-1.0, 'fredge%s'%l, 2)))
        for i,k,j,l in zip(fledge,bledge,bredge,fredge):
            mdb.models[modelName].Equation(name='E-3-fledge-bledge%s'%i, terms=((1.0, 'fledge%s'%i, 3), (-1.0, 'bledge%s'%k, 3)))
            mdb.models[modelName].Equation(name='E-3-bledge-bredge%s'%k, terms=((1.0, 'bledge%s'%k, 3), (-1.0, 'bredge%s'%j, 3),(-1.0, 'RP6', 3)))
            mdb.models[modelName].Equation(name='E-3-bredge-fredge%s'%j, terms=((1.0, 'bredge%s'%j, 3), (-1.0, 'fredge%s'%l, 3)))

        for i,k,j,l in zip(ltedge,lbedge,rbedge,rtedge):
            mdb.models[modelName].Equation(name='E-1-ltedge-lbedge%s'%i, terms=((1.0, 'ltedge%s'%i, 1), (-1.0, 'lbedge%s'%k, 1)))
            mdb.models[modelName].Equation(name='E-1-lbtedge-rbedge%s'%k, terms=((1.0, 'lbedge%s'%k, 1), (-1.0, 'rbedge%s'%j, 1)))
            mdb.models[modelName].Equation(name='E-1-rbedge-rtbedge%s'%j, terms=((1.0, 'rbedge%s'%j, 1), (-1.0, 'rtedge%s'%l, 1)))                                    
        for i,k,j,l in zip(ltedge,lbedge,rbedge,rtedge):
            mdb.models[modelName].Equation(name='E-2-ltedge-lbedge%s'%i, terms=((1.0, 'ltedge%s'%i, 2), (-1.0, 'lbedge%s'%k, 2),(-1.0, 'RP5', 2)))
            mdb.models[modelName].Equation(name='E-2-lbtedge-rbedge%s'%k, terms=((1.0, 'lbedge%s'%k, 2), (-1.0, 'rbedge%s'%j, 2)))
            mdb.models[modelName].Equation(name='E-2-rbedge-rtbedge%s'%j, terms=((1.0, 'rbedge%s'%j, 2), (-1.0, 'rtedge%s'%l, 2),(1.0, 'RP5', 2)))
        for i,k,j,l in zip(ltedge,lbedge,rbedge,rtedge):
            mdb.models[modelName].Equation(name='E-3-ltedge-lbedge%s'%i, terms=((1.0, 'ltedge%s'%i, 3), (-1.0, 'lbedge%s'%k, 3)))
            mdb.models[modelName].Equation(name='E-3-lbtedge-rbedge%s'%k, terms=((1.0, 'lbedge%s'%k, 3), (-1.0, 'rbedge%s'%j, 3),(-1.0, 'RP6', 3)))
            mdb.models[modelName].Equation(name='E-3-rbedge-rtbedge%s'%j, terms=((1.0, 'rbedge%s'%j, 3), (-1.0, 'rtedge%s'%l, 3)))

        #endregion

        #region Elastic modulus E11 ##
        for i in mdb.models[modelName].loads.keys():
                del mdb.models[modelName].loads[i]
        for i in mdb.models[modelName].boundaryConditions.keys():
                del mdb.models[modelName].boundaryConditions[i]


        region = a.sets['RP4']
        mdb.models[modelName].DisplacementBC(name='E11-1', createStepName='Step-1', 
            region=region, u1=Dispx, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
            amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
            localCsys=None)

        regionDef=mdb.models[modelName].rootAssembly.sets['c1']
        mdb.models[modelName].HistoryOutputRequest(name='H-Output-2', 
            createStepName='Step-1', variables=('RT', ), region=regionDef, 
            sectionPoints=DEFAULT, rebar=EXCLUDE)

        import os, glob

        mdb.Job(name='job-E11', model= modelName, description='', type=ANALYSIS, 
            atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
            memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
            explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
            modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
            scratch='', multiprocessingMode=DEFAULT, numCpus=CPUs, numDomains=CPUs, numGPUs=0)
        mdb.jobs['job-E11'].submit(consistencyChecking=OFF)
        
        mdb.jobs['job-E11'].waitForCompletion()
        o3 = session.openOdb(name='%s' % (path+'\job-E11.odb'))
        
        odb = session.odbs['%s' % (path+'\job-E11.odb')]

        session.viewports['Viewport: 1'].setValues(displayedObject=o3)
        odbName=session.viewports[session.currentViewportName].odbDisplay.name



        for i in session.xyDataObjects.keys():
            del session.xyDataObjects['%s' % (i)]

        session.odbData[odbName].setValues(activeFrames=(('Step-1', (1, )), ))
        session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('RF', 
            NODAL, ((COMPONENT, 'RF1'), )), ), nodeSets=('RP4', ))

        forceE11 = 0
        for i in session.xyDataObjects.keys():
            forceE11=forceE11+(session.xyDataObjects[i][0][1])

        stressE11 = abs(forceE11/(H*W))


        E11 = stressE11/(Dispx/L)                                


        for i in session.xyDataObjects.keys():
            del session.xyDataObjects['%s' % (i)]
        
        session.odbData[odbName].setValues(activeFrames=(('Step-1', (1, )), ))
        session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('U', 
            NODAL, ((COMPONENT, 'U1'), )), ), nodeSets=('C1','C2', ))
        
        C1U1new = session.xyDataObjects['U:U1 PI: %s N: %s' % (upperName,c1[0])][0][1] + coc1[(c1[0])][0]
        
        C2U1new = session.xyDataObjects['U:U1 PI: %s N: %s' % (upperName,c2[0])][0][1] + coc2[(c2[0])][0]
        Dis = abs(C1U1new - C2U1new)

        E11U1= abs(L - Dis)


        for i in session.xyDataObjects.keys():
            del session.xyDataObjects['%s' % (i)]
            

        session.odbData[odbName].setValues(activeFrames=(('Step-1', (1, )), ))
        session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('U', 
            NODAL, ((COMPONENT, 'U2'), )), ), nodeSets=('C1','C5', ))
        
        C1U2new = session.xyDataObjects['U:U2 PI: %s N: %s' % (upperName,c1[0])][0][1] + coc1[(c1[0])][1]
        
        C5U2new = session.xyDataObjects['U:U2 PI: %s N: %s' % (upperName,c5[0])][0][1] + coc5[(c5[0])][1]
        Dis = abs(C1U2new - C5U2new)

        E11U2= abs(H - Dis)

        for i in session.xyDataObjects.keys():
            del session.xyDataObjects['%s' % (i)]


        session.odbData[odbName].setValues(activeFrames=(('Step-1', (1, )), ))
        session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('U', 
            NODAL, ((COMPONENT, 'U3'), )), ), nodeSets=('C1','C4', ))
        
        C1U3new = session.xyDataObjects['U:U3 PI: %s N: %s' % (upperName,c1[0])][0][1] + coc1[(c1[0])][2]
        
        C4U3new = session.xyDataObjects['U:U3 PI: %s N: %s' % (upperName,c4[0])][0][1] + coc4[(c4[0])][2]
        Dis = abs(C1U3new - C4U3new)

        E11U3= abs(W - Dis)


        V12=(E11U2/H)/(E11U1/L)
        V13=(E11U3/W)/(E11U1/L)

        #endregion

        #region Elastic modulus E22 ##                                
        for i in mdb.models[modelName].loads.keys():
                del mdb.models[modelName].loads[i]
        for i in mdb.models[modelName].boundaryConditions.keys():
                del mdb.models[modelName].boundaryConditions[i]

        region = a.sets['RP5']
        mdb.models[modelName].DisplacementBC(name='E22-1', createStepName='Step-1', 
            region=region, u1=UNSET, u2=Dispy, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
            amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
            localCsys=None)


        regionDef=mdb.models[modelName].rootAssembly.sets['c1']
        mdb.models[modelName].HistoryOutputRequest(name='H-Output-2', 
            createStepName='Step-1', variables=('RT', ), region=regionDef, 
            sectionPoints=DEFAULT, rebar=EXCLUDE)

        import os, glob

        mdb.Job(name='job-E22', model= modelName, description='', type=ANALYSIS, 
            atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
            memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
            explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
            modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
            scratch='', multiprocessingMode=DEFAULT, numCpus=CPUs, numDomains=CPUs, numGPUs=0)
        mdb.jobs['job-E22'].submit(consistencyChecking=OFF)
        mdb.jobs['job-E22'].waitForCompletion()
        o3 = session.openOdb(name='%s' % (path+'\job-E22.odb'))
        
        odb = session.odbs['%s' % (path+'\job-E22.odb')]

        session.viewports['Viewport: 1'].setValues(displayedObject=o3)
        odbName=session.viewports[session.currentViewportName].odbDisplay.name


        for i in session.xyDataObjects.keys():
            del session.xyDataObjects['%s' % (i)]

        session.odbData[odbName].setValues(activeFrames=(('Step-1', (1, )), ))
        session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('RF', 
            NODAL, ((COMPONENT, 'RF2'), )), ), nodeSets=('RP5', ))

        forceE22 = 0
        for i in session.xyDataObjects.keys():
            forceE22=forceE22+(session.xyDataObjects[i][0][1])

        stressE22 = abs(forceE22/(W*L))


        E22 = stressE22/(Dispy/H)                                




        for i in session.xyDataObjects.keys():
            del session.xyDataObjects['%s' % (i)]


        
        session.odbData[odbName].setValues(activeFrames=(('Step-1', (1, )), ))
        session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('U', 
            NODAL, ((COMPONENT, 'U1'), )), ), nodeSets=('C1','C2', ))
        
        C1U1new = session.xyDataObjects['U:U1 PI: %s N: %s' % (upperName,c1[0])][0][1] + coc1[(c1[0])][0]
        
        C2U1new = session.xyDataObjects['U:U1 PI: %s N: %s' % (upperName,c2[0])][0][1] + coc2[(c2[0])][0]
        Dis = abs(C1U1new - C2U1new)

        E22U1= abs(L - Dis)


        for i in session.xyDataObjects.keys():
            del session.xyDataObjects['%s' % (i)]
            

        session.odbData[odbName].setValues(activeFrames=(('Step-1', (1, )), ))
        session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('U', 
            NODAL, ((COMPONENT, 'U2'), )), ), nodeSets=('C1','C5', ))
        
        C1U2new = session.xyDataObjects['U:U2 PI: %s N: %s' % (upperName,c1[0])][0][1] + coc1[(c1[0])][1]
        
        C5U2new = session.xyDataObjects['U:U2 PI: %s N: %s' % (upperName,c5[0])][0][1] + coc5[(c5[0])][1]
        Dis = abs(C1U2new - C5U2new)

        E22U2= abs(H - Dis)

        for i in session.xyDataObjects.keys():
            del session.xyDataObjects['%s' % (i)]


        session.odbData[odbName].setValues(activeFrames=(('Step-1', (1, )), ))
        session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('U', 
            NODAL, ((COMPONENT, 'U3'), )), ), nodeSets=('C1','C4', ))
        
        C1U3new = session.xyDataObjects['U:U3 PI: %s N: %s' % (upperName,c1[0])][0][1] + coc1[(c1[0])][2]
        
        C4U3new = session.xyDataObjects['U:U3 PI: %s N: %s' % (upperName,c4[0])][0][1] + coc4[(c4[0])][2]
        Dis = abs(C1U3new - C4U3new)

        E22U3= abs(W - Dis)


        V21=(E22U1/L)/(E22U2/H)
        V23=(E22U3/W)/(E22U2/H)
        #endregion

        #region Elastic modulus E33 ##
        if E33==True:
            for i in mdb.models[modelName].loads.keys():
                    del mdb.models[modelName].loads[i]
            for i in mdb.models[modelName].boundaryConditions.keys():
                    del mdb.models[modelName].boundaryConditions[i]


            region = a.sets['RP6']
            mdb.models[modelName].DisplacementBC(name='E33-1', createStepName='Step-1', 
                region=region, u1=UNSET, u2=UNSET, u3=Dispz, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
                amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
                localCsys=None)

            regionDef=mdb.models[modelName].rootAssembly.sets['c1']
            mdb.models[modelName].HistoryOutputRequest(name='H-Output-2', 
                createStepName='Step-1', variables=('RT', ), region=regionDef, 
                sectionPoints=DEFAULT, rebar=EXCLUDE)

            import os, glob

            mdb.Job(name='job-E33', model= modelName, description='', type=ANALYSIS, 
                atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
                memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
                explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
                modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
                scratch='', multiprocessingMode=DEFAULT, numCpus=CPUs, numDomains=CPUs, numGPUs=0)
            mdb.jobs['job-E33'].submit(consistencyChecking=OFF)
            mdb.jobs['job-E33'].waitForCompletion()
            o3 = session.openOdb(name='%s' % (path+'\job-E33.odb'))
            
            odb = session.odbs['%s' % (path+'\job-E33.odb')]

            session.viewports['Viewport: 1'].setValues(displayedObject=o3)
            odbName=session.viewports[session.currentViewportName].odbDisplay.name

            for i in session.xyDataObjects.keys():
                del session.xyDataObjects['%s' % (i)]

            session.odbData[odbName].setValues(activeFrames=(('Step-1', (1, )), ))
            session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('RF', 
                NODAL, ((COMPONENT, 'RF3'), )), ), nodeSets=('RP6', ))

            forceE33 = 0
            for i in session.xyDataObjects.keys():
                forceE33=forceE33+(session.xyDataObjects[i][0][1])

            stressE33 = abs(forceE33/(H*L))


            E33 = stressE33/(Dispz/W)                                

            for i in session.xyDataObjects.keys():
                del session.xyDataObjects['%s' % (i)]


            
            session.odbData[odbName].setValues(activeFrames=(('Step-1', (1, )), ))
            session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('U', 
                NODAL, ((COMPONENT, 'U1'), )), ), nodeSets=('C1','C2', ))
            
            C1U1new = session.xyDataObjects['U:U1 PI: %s N: %s' % (upperName,c1[0])][0][1] + coc1[(c1[0])][0]
            
            C2U1new = session.xyDataObjects['U:U1 PI: %s N: %s' % (upperName,c2[0])][0][1] + coc2[(c2[0])][0]
            Dis = abs(C1U1new - C2U1new)

            E33U1= abs(L - Dis)


            for i in session.xyDataObjects.keys():
                del session.xyDataObjects['%s' % (i)]
                

            session.odbData[odbName].setValues(activeFrames=(('Step-1', (1, )), ))
            session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('U', 
                NODAL, ((COMPONENT, 'U2'), )), ), nodeSets=('C1','C5', ))
            
            C1U2new = session.xyDataObjects['U:U2 PI: %s N: %s' % (upperName,c1[0])][0][1] + coc1[(c1[0])][1]
            
            C5U2new = session.xyDataObjects['U:U2 PI: %s N: %s' % (upperName,c5[0])][0][1] + coc5[(c5[0])][1]
            Dis = abs(C1U2new - C5U2new)

            E33U2= abs(H - Dis)

            for i in session.xyDataObjects.keys():
                del session.xyDataObjects['%s' % (i)]


            session.odbData[odbName].setValues(activeFrames=(('Step-1', (1, )), ))
            session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('U', 
                NODAL, ((COMPONENT, 'U3'), )), ), nodeSets=('C1','C4', ))
            
            C1U3new = session.xyDataObjects['U:U3 PI: %s N: %s' % (upperName,c1[0])][0][1] + coc1[(c1[0])][2]
            
            C4U3new = session.xyDataObjects['U:U3 PI: %s N: %s' % (upperName,c4[0])][0][1] + coc4[(c4[0])][2]
            Dis = abs(C1U3new - C4U3new)

            E33U3= abs(W - Dis)


            V31=(E33U1/L)/(E33U3/W)
            V32=(E33U2/H)/(E33U3/W)

        if E33==False:
                E33='N/A'
                V31='N/A'
                V32='N/A'
        #endregion


        #region Creating constraints for shear moduli ##
        for i in mdb.models[modelName].constraints.keys():
                del mdb.models[modelName].constraints[i]

        for i,k in zip(tops,bots):
            mdb.models[modelName].Equation(name='G-1-tops-bots%s'%i, terms=((1.0, 'tops%s'%i, 1), (-1.0, 'bots%s'%k, 1),(-1.0, 'RP4', 1)))
        for i,k in zip(tops,bots):
            mdb.models[modelName].Equation(name='G-2-tops-bots%s'%i, terms=((1.0, 'tops%s'%i, 2), (-1.0, 'bots%s'%k, 2),(-1.0, 'RP1', 2)))
        for i,k in zip(tops,bots):
            mdb.models[modelName].Equation(name='G-3-tops-bots%s'%i, terms=((1.0, 'tops%s'%i, 3), (-1.0, 'bots%s'%k, 3),(-1.0, 'RP6', 3)))

        for i,k in zip(lefts,rights):
            mdb.models[modelName].Equation(name='G-1-lefts-rights%s'%i, terms=((1.0, 'lefts%s'%i, 1), (-1.0, 'rights%s'%k, 1),(-1.0, 'RP5', 1)))
        for i,k in zip(lefts,rights):
            mdb.models[modelName].Equation(name='G-2-lefts-rights%s'%i, terms=((1.0, 'lefts%s'%i, 2), (-1.0, 'rights%s'%k, 2),(-1.0, 'RP6', 2)))
        for i,k in zip(lefts,rights):
            mdb.models[modelName].Equation(name='G-3-lefts-rights%s'%i, terms=((1.0, 'lefts%s'%i, 3), (-1.0, 'rights%s'%k, 3),(-1.0, 'RP2', 3)))

        for i,k in zip(fronts,backs):
            mdb.models[modelName].Equation(name='G-1-fronts-backs%s'%i, terms=((1.0, 'fronts%s'%i, 1), (-1.0, 'backs%s'%k, 1),(-1.0, 'RP3', 1)))
        for i,k in zip(fronts,backs):
            mdb.models[modelName].Equation(name='G-2-fronts-backs%s'%i, terms=((1.0, 'fronts%s'%i, 2), (-1.0, 'backs%s'%k, 2),(-1.0, 'RP4', 2)))
        for i,k in zip(fronts,backs):
            mdb.models[modelName].Equation(name='G-3-fronts-backs%s'%i, terms=((1.0, 'fronts%s'%i, 3), (-1.0, 'backs%s'%k, 3),(-1.0, 'RP5', 3)))


        mdb.models[modelName].Equation(name='G-1-c12', terms=((1.0, 'c6', 1), (-1.0, 'c2', 1),(1.0, 'RP4', 1)))
        mdb.models[modelName].Equation(name='G-1-c23', terms=((1.0, 'c2', 1), (-1.0, 'c3', 1),(-1.0, 'RP5', 1)))
        mdb.models[modelName].Equation(name='G-1-c34', terms=((1.0, 'c3', 1), (-1.0, 'c4', 1),(1.0, 'RP3', 1)))
        mdb.models[modelName].Equation(name='G-1-c45', terms=((1.0, 'c4', 1), (-1.0, 'c8', 1),(-1.0, 'RP4', 1)))
        mdb.models[modelName].Equation(name='G-1-c56', terms=((1.0, 'c8', 1), (-1.0, 'c5', 1),(1.0, 'RP5', 1)))
        mdb.models[modelName].Equation(name='G-1-c67', terms=((1.0, 'c5', 1), (-1.0, 'c1', 1),(1.0, 'RP4', 1)))
        mdb.models[modelName].Equation(name='G-1-c78', terms=((1.0, 'c1', 1), (-1.0, 'c7', 1),(-1.0, 'RP3', 1),(-1.0, 'RP4', 1),(-1.0, 'RP5', 1)))

            
        mdb.models[modelName].Equation(name='G-2-c12', terms=((1.0, 'c6', 2), (-1.0, 'c2', 2),(1.0, 'RP1', 2)))
        mdb.models[modelName].Equation(name='G-2-c23', terms=((1.0, 'c2', 2), (-1.0, 'c3', 2),(-1.0, 'RP6', 2)))
        mdb.models[modelName].Equation(name='G-2-c34', terms=((1.0, 'c3', 2), (-1.0, 'c4', 2),(1.0, 'RP4', 2)))
        mdb.models[modelName].Equation(name='G-2-c45', terms=((1.0, 'c4', 2), (-1.0, 'c8', 2),(-1.0, 'RP1', 2)))
        mdb.models[modelName].Equation(name='G-2-c56', terms=((1.0, 'c8', 2), (-1.0, 'c5', 2),(1.0, 'RP6', 2)))
        mdb.models[modelName].Equation(name='G-2-c67', terms=((1.0, 'c5', 2), (-1.0, 'c1', 2),(1.0, 'RP1', 2)))
        mdb.models[modelName].Equation(name='G-2-c78', terms=((1.0, 'c1', 2), (-1.0, 'c7', 2),(-1.0, 'RP1', 2),(-1.0, 'RP4', 2),(-1.0, 'RP6', 2)))


        mdb.models[modelName].Equation(name='G-3-c12', terms=((1.0, 'c6', 3), (-1.0, 'c2', 3),(1.0, 'RP6', 3)))
        mdb.models[modelName].Equation(name='G-3-c23', terms=((1.0, 'c2', 3), (-1.0, 'c3', 3),(-1.0, 'RP2', 3)))
        mdb.models[modelName].Equation(name='G-3-c34', terms=((1.0, 'c3', 3), (-1.0, 'c4', 3),(1.0, 'RP5', 3)))
        mdb.models[modelName].Equation(name='G-3-c45', terms=((1.0, 'c4', 3), (-1.0, 'c8', 3),(-1.0, 'RP6', 3)))
        mdb.models[modelName].Equation(name='G-3-c56', terms=((1.0, 'c8', 3), (-1.0, 'c5', 3),(1.0, 'RP2', 3)))
        mdb.models[modelName].Equation(name='G-3-c67', terms=((1.0, 'c5', 3), (-1.0, 'c1', 3),(1.0, 'RP6', 3)))
        mdb.models[modelName].Equation(name='G-3-c78', terms=((1.0, 'c1', 3), (-1.0, 'c7', 3),(-1.0, 'RP2', 3),(-1.0, 'RP5', 3),(-1.0, 'RP6', 3)))


        for i,k,j,l in zip(ftedge,btedge,bbedge,fbedge):
            mdb.models[modelName].Equation(name='G-1-ftedge-btedge%s'%i, terms=((1.0, 'ftedge%s'%i, 1), (-1.0, 'btedge%s'%k, 1),(-1.0, 'RP3', 1)))
            mdb.models[modelName].Equation(name='G-1-btedge-bbedge%s'%k, terms=((1.0, 'btedge%s'%k, 1), (-1.0, 'bbedge%s'%j, 1),(-1.0, 'RP4', 1)))
            mdb.models[modelName].Equation(name='G-1-bbedge-fbedge%s'%j, terms=((1.0, 'bbedge%s'%j, 1), (-1.0, 'fbedge%s'%l, 1),(1.0, 'RP3', 1)))
        for i,k,j,l in zip(ftedge,btedge,bbedge,fbedge):
            mdb.models[modelName].Equation(name='G-2-ftedge-btedge%s'%i, terms=((1.0, 'ftedge%s'%i, 2), (-1.0, 'btedge%s'%k, 2),(-1.0, 'RP4', 2)))
            mdb.models[modelName].Equation(name='G-2-btedge-bbedge%s'%k, terms=((1.0, 'btedge%s'%k, 2), (-1.0, 'bbedge%s'%j, 2),(-1.0, 'RP1', 2)))
            mdb.models[modelName].Equation(name='G-2-bbedge-fbedge%s'%j, terms=((1.0, 'bbedge%s'%j, 2), (-1.0, 'fbedge%s'%l, 2),(1.0, 'RP4', 2)))
        for i,k,j,l in zip(ftedge,btedge,bbedge,fbedge):
            mdb.models[modelName].Equation(name='G-3-ftedge-btedge%s'%i, terms=((1.0, 'ftedge%s'%i, 3), (-1.0, 'btedge%s'%k, 3),(-1.0, 'RP5', 3)))
            mdb.models[modelName].Equation(name='G-3-btedge-bbedge%s'%k, terms=((1.0, 'btedge%s'%k, 3), (-1.0, 'bbedge%s'%j, 3),(-1.0, 'RP6', 3)))
            mdb.models[modelName].Equation(name='G-3-bbedge-fbedge%s'%j, terms=((1.0, 'bbedge%s'%j, 3), (-1.0, 'fbedge%s'%l, 3),(1.0, 'RP5', 3)))


        for i,k,j,l in zip(fledge,bledge,bredge,fredge):
            mdb.models[modelName].Equation(name='G-1-fledge-bledge%s'%i, terms=((1.0, 'fledge%s'%i, 1), (-1.0, 'bledge%s'%k, 1),(-1.0, 'RP3', 1)))
            mdb.models[modelName].Equation(name='G-1-bledge-bredge%s'%k, terms=((1.0, 'bledge%s'%k, 1), (-1.0, 'bredge%s'%j, 1),(-1.0, 'RP5', 1)))
            mdb.models[modelName].Equation(name='G-1-bredge-fredge%s'%j, terms=((1.0, 'bredge%s'%j, 1), (-1.0, 'fredge%s'%l, 1),(1.0, 'RP3', 1)))
        for i,k,j,l in zip(fledge,bledge,bredge,fredge):
            mdb.models[modelName].Equation(name='G-2-fledge-bledge%s'%i, terms=((1.0, 'fledge%s'%i, 2), (-1.0, 'bledge%s'%k, 2),(-1.0, 'RP4', 2)))
            mdb.models[modelName].Equation(name='G-2-bledge-bredge%s'%k, terms=((1.0, 'bledge%s'%k, 2), (-1.0, 'bredge%s'%j, 2),(-1.0, 'RP6', 2)))
            mdb.models[modelName].Equation(name='G-2-bredge-fredge%s'%j, terms=((1.0, 'bredge%s'%j, 2), (-1.0, 'fredge%s'%l, 2),(1.0, 'RP4', 2)))
        for i,k,j,l in zip(fledge,bledge,bredge,fredge):
            mdb.models[modelName].Equation(name='G-3-fledge-bledge%s'%i, terms=((1.0, 'fledge%s'%i, 3), (-1.0, 'bledge%s'%k, 3),(-1.0, 'RP5', 3)))
            mdb.models[modelName].Equation(name='G-3-bledge-bredge%s'%k, terms=((1.0, 'bledge%s'%k, 3), (-1.0, 'bredge%s'%j, 3),(-1.0, 'RP2', 3)))
            mdb.models[modelName].Equation(name='G-3-bredge-fredge%s'%j, terms=((1.0, 'bredge%s'%j, 3), (-1.0, 'fredge%s'%l, 3),(1.0, 'RP5', 3)))
            

        for i,k,j,l in zip(ltedge,lbedge,rbedge,rtedge):
            mdb.models[modelName].Equation(name='G-1-ltedge-lbedge%s'%i, terms=((1.0, 'ltedge%s'%i, 1), (-1.0, 'lbedge%s'%k, 1),(-1.0, 'RP4', 1)))
            mdb.models[modelName].Equation(name='G-1-lbtedge-rbedge%s'%k, terms=((1.0, 'lbedge%s'%k, 1), (-1.0, 'rbedge%s'%j, 1),(-1.0, 'RP5', 1)))
            mdb.models[modelName].Equation(name='G-1-rbedge-rtbedge%s'%j, terms=((1.0, 'rbedge%s'%j, 1), (-1.0, 'rtedge%s'%l, 1),(1.0, 'RP4', 1)))                                    
        for i,k,j,l in zip(ltedge,lbedge,rbedge,rtedge):
            mdb.models[modelName].Equation(name='G-2-ltedge-lbedge%s'%i, terms=((1.0, 'ltedge%s'%i, 2), (-1.0, 'lbedge%s'%k, 2),(-1.0, 'RP1', 2)))
            mdb.models[modelName].Equation(name='G-2-lbtedge-rbedge%s'%k, terms=((1.0, 'lbedge%s'%k, 2), (-1.0, 'rbedge%s'%j, 2),(-1.0, 'RP6', 2)))
            mdb.models[modelName].Equation(name='G-2-rbedge-rtbedge%s'%j, terms=((1.0, 'rbedge%s'%j, 2), (-1.0, 'rtedge%s'%l, 2),(1.0, 'RP1', 2)))
        for i,k,j,l in zip(ltedge,lbedge,rbedge,rtedge):
            mdb.models[modelName].Equation(name='G-3-ltedge-lbedge%s'%i, terms=((1.0, 'ltedge%s'%i, 3), (-1.0, 'lbedge%s'%k, 3),(-1.0, 'RP6', 3)))
            mdb.models[modelName].Equation(name='G-3-lbtedge-rbedge%s'%k, terms=((1.0, 'lbedge%s'%k, 3), (-1.0, 'rbedge%s'%j, 3),(-1.0, 'RP2', 3)))
            mdb.models[modelName].Equation(name='G-3-rbedge-rtbedge%s'%j, terms=((1.0, 'rbedge%s'%j, 3), (-1.0, 'rtedge%s'%l, 3),(1.0, 'RP6', 3)))
        #endregion

        #region Shear modulus G12 ##
        for i in mdb.models[modelName].loads.keys():
                del mdb.models[modelName].loads[i]
        for i in mdb.models[modelName].boundaryConditions.keys():
                del mdb.models[modelName].boundaryConditions[i]

        region = a.sets['RP4']
        mdb.models[modelName].DisplacementBC(name='G12-1', createStepName='Step-1', 
            region=region, u1=Dispx, u2=Dispy, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
            amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
            localCsys=None)


        region = a.sets['RP5']
        mdb.models[modelName].DisplacementBC(name='G12-2', createStepName='Step-1', 
            region=region, u1=0, u2=0, u3=0, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
            amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
            localCsys=None)

        region = a.sets['RP6']
        mdb.models[modelName].DisplacementBC(name='G12-3', createStepName='Step-1', 
            region=region, u1=0, u2=0, u3=0, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
            amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
            localCsys=None)

        regionDef=mdb.models[modelName].rootAssembly.sets['c1']
        mdb.models[modelName].HistoryOutputRequest(name='H-Output-2', 
            createStepName='Step-1', variables=('RT', ), region=regionDef, 
            sectionPoints=DEFAULT, rebar=EXCLUDE)

        import os, glob

        mdb.Job(name='job-G12', model= modelName, description='', type=ANALYSIS, 
            atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
            memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
            explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
            modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
            scratch='', multiprocessingMode=DEFAULT, numCpus=CPUs, numDomains=CPUs, numGPUs=0)
        mdb.jobs['job-G12'].submit(consistencyChecking=OFF)

        mdb.jobs['job-G12'].waitForCompletion()


        o3 = session.openOdb(name='%s' % (path+'\job-G12.odb'))


        odb = session.odbs['%s' % (path+'\job-G12.odb')]

        session.viewports['Viewport: 1'].setValues(displayedObject=o3)
        odbName=session.viewports[session.currentViewportName].odbDisplay.name


        for i in session.xyDataObjects.keys():
            del session.xyDataObjects['%s' % (i)]

        session.odbData[odbName].setValues(activeFrames=(('Step-1', (1, )), ))
        session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('RF', 
            NODAL, ((COMPONENT, 'RF1'), )), ), nodeSets=('RP4', ))

        forceG12 = 0
        for i in session.xyDataObjects.keys():
            forceG12=forceG12+(session.xyDataObjects[i][0][1])

        stressG12 = abs(forceG12/(L*W))

        G12 = stressG12/((Dispx/H)+(Dispy/L))
        #endregion

        #region Shear modulus G13 ##                                
        for i in mdb.models[modelName].loads.keys():
                del mdb.models[modelName].loads[i]
        for i in mdb.models[modelName].boundaryConditions.keys():
                del mdb.models[modelName].boundaryConditions[i]


        region = a.sets['RP5']
        mdb.models[modelName].DisplacementBC(name='G13-1', createStepName='Step-1', 
            region=region, u1=Dispx, u2=UNSET, u3=Dispz, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
            amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
            localCsys=None)

        region = a.sets['RP4']
        mdb.models[modelName].DisplacementBC(name='G13-2', createStepName='Step-1', 
            region=region, u1=0, u2=0, u3=0, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
            amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
            localCsys=None)

        region = a.sets['RP6']
        mdb.models[modelName].DisplacementBC(name='G13-3', createStepName='Step-1', 
            region=region, u1=0, u2=0, u3=0, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
            amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
            localCsys=None)


        regionDef=mdb.models[modelName].rootAssembly.sets['c1']
        mdb.models[modelName].HistoryOutputRequest(name='H-Output-2', 
            createStepName='Step-1', variables=('RT', ), region=regionDef, 
            sectionPoints=DEFAULT, rebar=EXCLUDE)

        import os, glob

        mdb.Job(name='job-G13', model= modelName, description='', type=ANALYSIS, 
            atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
            memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
            explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
            modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
            scratch='', multiprocessingMode=DEFAULT, numCpus=CPUs, numDomains=CPUs, numGPUs=0)
        mdb.jobs['job-G13'].submit(consistencyChecking=OFF)

        mdb.jobs['job-G13'].waitForCompletion()


        o3 = session.openOdb(name='%s' % (path+'\job-G13.odb'))


        odb = session.odbs['%s' % (path+'\job-G13.odb')]

        session.viewports['Viewport: 1'].setValues(displayedObject=o3)
        odbName=session.viewports[session.currentViewportName].odbDisplay.name



        for i in session.xyDataObjects.keys():
            del session.xyDataObjects['%s' % (i)]

        session.odbData[odbName].setValues(activeFrames=(('Step-1', (1, )), ))
        session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('RF', 
            NODAL, ((COMPONENT, 'RF1'), )), ), nodeSets=('RP5', ))

        forceG13 = 0
        for i in session.xyDataObjects.keys():
            forceG13=forceG13+(session.xyDataObjects[i][0][1])

        stressG13 = abs(forceG13/(H*L))


        G13 = stressG13/((Dispx/W)+(Dispz/L))
        #endregion

        #region Shear modulus G23 ##
        for i in mdb.models[modelName].loads.keys():
                del mdb.models[modelName].loads[i]
        for i in mdb.models[modelName].boundaryConditions.keys():
                del mdb.models[modelName].boundaryConditions[i]


        region = a.sets['RP6']
        mdb.models[modelName].DisplacementBC(name='G23-1', createStepName='Step-1', 
            region=region, u1=UNSET, u2=Dispy, u3=Dispz, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
            amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
            localCsys=None)

        region = a.sets['RP4']
        mdb.models[modelName].DisplacementBC(name='G23-2', createStepName='Step-1', 
            region=region, u1=0, u2=0, u3=0, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
            amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
            localCsys=None)

        region = a.sets['RP5']
        mdb.models[modelName].DisplacementBC(name='G23-3', createStepName='Step-1', 
            region=region, u1=0, u2=0, u3=0, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
            amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
            localCsys=None)


        regionDef=mdb.models[modelName].rootAssembly.sets['c1']
        mdb.models[modelName].HistoryOutputRequest(name='H-Output-2', 
            createStepName='Step-1', variables=('RT', ), region=regionDef, 
            sectionPoints=DEFAULT, rebar=EXCLUDE)

        import os, glob


        mdb.Job(name='job-G23', model= modelName, description='', type=ANALYSIS, 
            atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
            memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
            explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
            modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
            scratch='', multiprocessingMode=DEFAULT, numCpus=CPUs, numDomains=CPUs, numGPUs=0)
        mdb.jobs['job-G23'].submit(consistencyChecking=OFF)

        mdb.jobs['job-G23'].waitForCompletion()

        o3 = session.openOdb(name='%s' % (path+'\job-G23.odb'))

        odb = session.odbs['%s' % (path+'\job-G23.odb')]

        session.viewports['Viewport: 1'].setValues(displayedObject=o3)
        odbName=session.viewports[session.currentViewportName].odbDisplay.name

        for i in session.xyDataObjects.keys():
            del session.xyDataObjects['%s' % (i)]

        session.odbData[odbName].setValues(activeFrames=(('Step-1', (1, )), ))
        session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('RF', 
            NODAL, ((COMPONENT, 'RF2'), )), ), nodeSets=('RP6', ))

        forceG23 = 0
        for i in session.xyDataObjects.keys():
            forceG23=forceG23+(session.xyDataObjects[i][0][1])

        stressG23 = abs(forceG23/(L*H))


        G23 = stressG23/((Dispy/W)+(Dispz/H))

        #endregion


    #################################################
        print ('----------------------------------------------------')
        print ('The homogenised elastic properties:')
        print ('E11=%s Stress units' % (E11))
        print ('V12=%s ratio' % (V12))
        print ('V13=%s ratio' % (V13))
        print ('E22=%s Stress units' % (E22))
        print ('V21=%s ratio' % (V21))
        print ('V23=%s ratio' % (V23))
        print ('E33=%s Stress units' % (E33))
        print ('V31=%s ratio' % (V31))
        print ('V32=%s ratio' % (V32))
        print ('G12=%s Stress units' % (G12))
        print ('G13=%s Stress units' % (G13))
        print ('G23=%s Stress units' % (G23))
        print ('----------------------------------------------------')
        print ('Processing duration %s seconds' % (time.time()-start))
        print ('----------------------------------------------------')
        
        filename = ('%s_elastic_properties.txt' % part)
        print ('The homogenised elastic properties are saved in ABAQUS Work Directory under %s' % filename)
        f = open(filename,'w')
        f.write('{0:^10}{1:^20}{2:^20}\n'.format('Property','Value','Unit'))
        f.write('{0:^10}{1:^20}{2:^20}\n'.format('E11',E11,'Stress units'))
        f.write('{0:^10}{1:^20}{2:^20}\n'.format('V12',V12,'ratio'))
        f.write('{0:^10}{1:^20}{2:^20}\n'.format('V13',V13,'ratio'))
        f.write('{0:^10}{1:^20}{2:^20}\n'.format('E22',E22,'Stress units'))
        f.write('{0:^10}{1:^20}{2:^20}\n'.format('V21',V21,'ratio'))
        f.write('{0:^10}{1:^20}{2:^20}\n'.format('V23',V23,'ratio'))
        f.write('{0:^10}{1:^20}{2:^20}\n'.format('E33',E33,'Stress units'))
        f.write('{0:^10}{1:^20}{2:^20}\n'.format('V31',V31,'ratio'))
        f.write('{0:^10}{1:^20}{2:^20}\n'.format('V32',V32,'ratio'))
        f.write('{0:^10}{1:^20}{2:^20}\n'.format('G12',G12,'Stress units'))
        f.write('{0:^10}{1:^20}{2:^20}\n'.format('G13',G13,'Stress units'))
        f.write('{0:^10}{1:^20}{2:^20}\n'.format('G23',G23,'Stress units'))               
        f.write ('processing duration %s Seconds' % (time.time()-start))

        f.close()

        filename = ('%s_elastic_properties(easycopy).txt' % part)
        f = open(filename,'w')
        f.write('{0:^10}\n'.format(E11))
        f.write('{0:^10}\n'.format(E22))
        f.write('{0:^10}\n'.format(E33))
        f.write('{0:^10}\n'.format(G12))
        f.write('{0:^10}\n'.format(G13))
        f.write('{0:^10}\n'.format(G23))
        f.write('{0:^10}\n'.format(V12))
        f.write('{0:^10}\n'.format(V13))
        f.write('{0:^10}\n'.format(V21))
        f.write('{0:^10}\n'.format(V23))
        f.write('{0:^10}\n'.format(V31))
        f.write('{0:^10}\n'.format(V32))                 
        f.write ('{0:^10}\n' .format((time.time()-start)))

        f.close()

        print ('----------------------------------------------------')
        for i in session.xyDataObjects.keys():
            del session.xyDataObjects['%s' % (i)]
        print ('--------- End of EasyPBC ---------')

        if len(session.odbData.keys()) >= 1:                              
                odb.close(odb, write=TRUE)
                a = mdb.models[modelName].rootAssembly
                session.viewports['Viewport: 1'].setValues(displayedObject=a)
    #endregion


    #============EasyPBC================
    a = mdb.models['Model-1'].rootAssembly
    a.regenerate()

    feasypbc(part='Model-1', inst='Part-1-1', meshsens=1E-07, CPU=CPUnum, E33=False)
    
    #=================================
    #endregion
    
    #region bending test
    #step
    mdb.models['Model-2'].StaticStep(name='Step-1', previous='Initial')


    #=============symmetry constraint method=============

    #===========job-BendingX===========
    #XsymmBC
    a = mdb.models['Model-2'].rootAssembly
    e1 = a.instances['Part-1-1'].edges
    edges1 = e1.getSequenceFromMask(mask=('[#8 #200 ]', ), )
    region = regionToolset.Region(edges=edges1)
    mdb.models['Model-2'].XsymmBC(name='BC-Xsymm', createStepName='Step-1', 
        region=region, localCsys=None)
    #YsymmBC
    a = mdb.models['Model-2'].rootAssembly
    e1 = a.instances['Part-1-1'].edges
    edges1 = e1.getSequenceFromMask(mask=('[#31249891 #129 ]', ), )
    v1 = a.instances['Part-1-1'].vertices
    verts1 = v1.getSequenceFromMask(mask=('[#154d54d3 ]', ), )
    region = regionToolset.Region(vertices=verts1, edges=edges1)
    mdb.models['Model-2'].YsymmBC(name='BC-Ysymm', createStepName='Step-1', 
        region=region, localCsys=None)
    #RP
    a = mdb.models['Model-2'].rootAssembly
    a.ReferencePoint(point=(length, width/2, h[j]/2))   

    #coupling-kinematic
    a = mdb.models['Model-2'].rootAssembly
    r1 = a.referencePoints
    refPoints1=(r1[6], )
    region1=a.Set(referencePoints=refPoints1, name='m_Set-rp1')
    a = mdb.models['Model-2'].rootAssembly
    e1 = a.instances['Part-1-1'].edges
    edges1 = e1.getSequenceFromMask(mask=('[#0 #82 ]', ), )
    region2=regionToolset.Region(edges=edges1)
    mdb.models['Model-2'].Coupling(name='Constraint-1', controlPoint=region1, 
        surface=region2, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, 
        localCsys=None, u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)
    #Load
    a = mdb.models['Model-2'].rootAssembly
    region = a.sets['m_Set-rp1']
    mdb.models['Model-2'].Moment(name='Load-1', createStepName='Step-1', 
        region=region, cm2=-width, distributionType=UNIFORM, field='', localCsys=None)    
        
    #job
    mdb.Job(name='Job-BendingX', model='Model-2', description='', type=ANALYSIS, 
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
        memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
        scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=CPUnum, 
        numDomains=CPUnum, numGPUs=0)

    mdb.jobs['Job-BendingX'].submit(consistencyChecking=OFF) 


    #===========job-BendingY===========
    #rp
    a = mdb.models['Model-2'].rootAssembly
    a.features['RP-1'].setValues(xValue=length/2, yValue=width)
    a = mdb.models['Model-2'].rootAssembly
    a.regenerate()

    #coupling
    a = mdb.models['Model-2'].rootAssembly
    e1 = a.instances['Part-1-1'].edges
    edges1 = e1.getSequenceFromMask(mask=('[#8c924644 #454 ]', ), )
    v1 = a.instances['Part-1-1'].vertices
    verts1 = v1.getSequenceFromMask(mask=('[#2ab2ab2c ]', ), )
    r2 = a.referencePoints
    refPoints2=(r2[6], )
    region2=regionToolset.Region(vertices=verts1, edges=edges1, 
        referencePoints=refPoints2)
    mdb.models['Model-2'].constraints['Constraint-1'].setValues(surface=region2)
    #Load
    mdb.models['Model-2'].loads['Load-1'].setValues(cm1=length,cm2=UNSET,
        distributionType=UNIFORM, field='')

    #job
    mdb.Job(name='Job-BendingY', model='Model-2', description='', type=ANALYSIS, 
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
        memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
        scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=CPUnum, 
        numDomains=CPUnum, numGPUs=0)

    mdb.jobs['Job-BendingY'].submit(consistencyChecking=OFF)  

    #===========job-BendingXY===========

    #rp
    a = mdb.models['Model-2'].rootAssembly
    del a.features['RP-1']

    a = mdb.models['Model-2'].rootAssembly
    e11 = a.instances['Part-1-1'].edges
    a.ReferencePoint(point=a.instances['Part-1-1'].InterestingPoint(edge=e11[39], 
        rule=MIDDLE))
    a = mdb.models['Model-2'].rootAssembly
    e1 = a.instances['Part-1-1'].edges
    a.ReferencePoint(point=a.instances['Part-1-1'].InterestingPoint(edge=e1[33], 
        rule=MIDDLE))
    a = mdb.models['Model-2'].rootAssembly
    v11 = a.instances['Part-1-1'].vertices
    a.ReferencePoint(point=v11[13])
    a = mdb.models['Model-2'].rootAssembly
    e11 = a.instances['Part-1-1'].edges
    a.ReferencePoint(point=a.instances['Part-1-1'].InterestingPoint(edge=e11[36], 
        rule=MIDDLE))


    #coupling
    del mdb.models['Model-2'].constraints['Constraint-1']
    del mdb.models['Model-2'].rootAssembly.sets['m_Set-rp1']

    a = mdb.models['Model-2'].rootAssembly
    r1 = a.referencePoints
    refPoints1=(r1[10], )
    region1=a.Set(referencePoints=refPoints1, name='m_Set-rp1')
    a = mdb.models['Model-2'].rootAssembly
    e1 = a.instances['Part-1-1'].edges
    edges1 = e1.getSequenceFromMask(mask=('[#0 #80 ]', ), )
    region2=regionToolset.Region(edges=edges1)
    mdb.models['Model-2'].Coupling(name='Constraint-1', controlPoint=region1, 
        surface=region2, influenceRadius=WHOLE_SURFACE, couplingType=DISTRIBUTING, 
        weightingMethod=UNIFORM, localCsys=None, u1=ON, u2=ON, u3=ON, ur1=ON, 
        ur2=ON, ur3=ON)
        
    a = mdb.models['Model-2'].rootAssembly
    r1 = a.referencePoints
    refPoints1=(r1[11], )
    region1=a.Set(referencePoints=refPoints1, name='m_Set-rp2')
    a = mdb.models['Model-2'].rootAssembly
    e1 = a.instances['Part-1-1'].edges
    edges1 = e1.getSequenceFromMask(mask=('[#0 #2 ]', ), )
    region2=regionToolset.Region(edges=edges1)
    mdb.models['Model-2'].Coupling(name='Constraint-2', controlPoint=region1, 
        surface=region2, influenceRadius=WHOLE_SURFACE, couplingType=DISTRIBUTING, 
        weightingMethod=UNIFORM, localCsys=None, u1=ON, u2=ON, u3=ON, ur1=ON, 
        ur2=ON, ur3=ON)    
    
    a = mdb.models['Model-2'].rootAssembly
    r1 = a.referencePoints
    refPoints1=(r1[12], )
    region1=a.Set(referencePoints=refPoints1, name='m_Set-rp3')
    a = mdb.models['Model-2'].rootAssembly
    e1 = a.instances['Part-1-1'].edges
    edges1 = e1.getSequenceFromMask(mask=('[#0 #440 ]', ), )
    region2=regionToolset.Region(edges=edges1)
    mdb.models['Model-2'].Coupling(name='Constraint-3', controlPoint=region1, 
        surface=region2, influenceRadius=WHOLE_SURFACE, couplingType=DISTRIBUTING, 
        weightingMethod=UNIFORM, localCsys=None, u1=ON, u2=ON, u3=ON, ur1=ON, 
        ur2=ON, ur3=ON)
    

        
    a = mdb.models['Model-2'].rootAssembly
    r1 = a.referencePoints
    refPoints1=(r1[13], )
    region1=a.Set(referencePoints=refPoints1, name='m_Set-rp4')
    a = mdb.models['Model-2'].rootAssembly
    e1 = a.instances['Part-1-1'].edges
    edges1 = e1.getSequenceFromMask(mask=('[#0 #10 ]', ), )
    region2=regionToolset.Region(edges=edges1)
    mdb.models['Model-2'].Coupling(name='Constraint-4', controlPoint=region1, 
        surface=region2, influenceRadius=WHOLE_SURFACE, couplingType=DISTRIBUTING, 
        weightingMethod=UNIFORM, localCsys=None, u1=ON, u2=ON, u3=ON, ur1=ON, 
        ur2=ON, ur3=ON)



    #BC
    mdb.models['Model-2'].boundaryConditions['BC-Xsymm'].setValues(typeName=XASYMM)
    mdb.models['Model-2'].boundaryConditions['BC-Ysymm'].setValues(typeName=YASYMM)

    a = mdb.models['Model-2'].rootAssembly
    region = a.sets['m_Set-rp1']
    mdb.models['Model-2'].DisplacementBC(name='BC-3', createStepName='Step-1', 
        region=region, u1=0.0, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
        amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
        localCsys=None)
    a = mdb.models['Model-2'].rootAssembly
    region = a.sets['m_Set-rp2']
    mdb.models['Model-2'].DisplacementBC(name='BC-4', createStepName='Step-1', 
        region=region, u1=0.0, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
        amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
        localCsys=None)
    a = mdb.models['Model-2'].rootAssembly
    region = a.sets['m_Set-rp3']
    mdb.models['Model-2'].DisplacementBC(name='BC-5', createStepName='Step-1', 
        region=region, u1=UNSET, u2=0.0, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
        amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
        localCsys=None)
    a = mdb.models['Model-2'].rootAssembly
    region = a.sets['m_Set-rp4']
    mdb.models['Model-2'].DisplacementBC(name='BC-6', createStepName='Step-1', 
        region=region, u1=UNSET, u2=0.0, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
        amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
        localCsys=None)
    #Load
    del mdb.models['Model-2'].loads['Load-1']

    a = mdb.models['Model-2'].rootAssembly
    region = a.sets['m_Set-rp1']
    mdb.models['Model-2'].ConcentratedForce(name='Load-1', createStepName='Step-1', 
        region=region, cf2=-1*width/h[j], distributionType=UNIFORM, field='', localCsys=None)

    a = mdb.models['Model-2'].rootAssembly
    region = a.sets['m_Set-rp2']
    mdb.models['Model-2'].ConcentratedForce(name='Load-2', createStepName='Step-1', 
        region=region, cf2=1*width/h[j], distributionType=UNIFORM, field='', 
        localCsys=None)

    a = mdb.models['Model-2'].rootAssembly
    region = a.sets['m_Set-rp3']
    mdb.models['Model-2'].ConcentratedForce(name='Load-3', createStepName='Step-1', 
        region=region, cf1=-1*length/h[j], distributionType=UNIFORM, field='', localCsys=None)

    a = mdb.models['Model-2'].rootAssembly
    region = a.sets['m_Set-rp4']
    mdb.models['Model-2'].ConcentratedForce(name='Load-4', createStepName='Step-1', 
        region=region, cf1=1*length/h[j], distributionType=UNIFORM, field='', 
        localCsys=None)

    
    #job
    mdb.Job(name='Job-BendingXY', model='Model-2', description='', type=ANALYSIS, 
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
        memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
        scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=CPUnum, 
        numDomains=CPUnum, numGPUs=0)
    mdb.jobs['Job-BendingXY'].submit(consistencyChecking=OFF) 

    #endregion

    #region extract node displacement
    #======extract node displacement===========
    # waitForCompletion
    mdb.jobs['Job-BendingX'].waitForCompletion()
    mdb.jobs['Job-BendingY'].waitForCompletion()  
    mdb.jobs['Job-BendingXY'].waitForCompletion()  

    from odbAccess import openOdb
    BendingX_UR = {}
    BendingY_UR = {}
    BendingXY_UR = {}



    #BendingX
    odb_path1 = 'Job-BendingX.odb'  
    odb1 = openOdb(odb_path1)
    ur_field1 = odb1.steps['Step-1'].frames[-1].fieldOutputs['UR']
    for index, node_label in enumerate(interested_nodes, start=1):
        BendingX_UR[index] = {'Node Label': node_label, 'UR':[0, 0, 0]}
    for disp_value in ur_field1.values:
        node_label = disp_value.nodeLabel
        if node_label in interested_nodes:
            ur_1 = disp_value.data[0]
            ur_2 = disp_value.data[1]
            ur_3 = disp_value.data[2]
            node_index = interested_nodes.index(node_label) + 1
            BendingX_UR[node_index]['UR'] = [ur_1, ur_2, ur_3]
    odb1.close()


    #BendingY
    odb_path2 = 'Job-BendingY.odb'  
    odb2 = openOdb(odb_path2)
    ur_field2 = odb2.steps['Step-1'].frames[-1].fieldOutputs['UR']
    for index, node_label in enumerate(interested_nodes, start=1):
        BendingY_UR[index] = {'Node Label': node_label, 'UR':[0, 0, 0]}
    for disp_value in ur_field2.values:
        node_label = disp_value.nodeLabel
        if node_label in interested_nodes:
            ur_1 = disp_value.data[0]
            ur_2 = disp_value.data[1]
            ur_3 = disp_value.data[2]
            node_index = interested_nodes.index(node_label) + 1
            BendingY_UR[node_index]['UR'] = [ur_1, ur_2, ur_3]
    odb2.close()  

    #BendingXY
    odb_path3 = 'Job-BendingXY.odb'  
    odb3 = openOdb(odb_path3)
    ur_field3 = odb3.steps['Step-1'].frames[-1].fieldOutputs['UR']
    for index, node_label in enumerate(interested_nodes, start=1):
        BendingXY_UR[index] = {'Node Label': node_label, 'UR':[0, 0, 0]}
    for disp_value in ur_field3.values:
        node_label = disp_value.nodeLabel
        if node_label in interested_nodes:
            ur_1 = disp_value.data[0]
            ur_2 = disp_value.data[1]
            ur_3 = disp_value.data[2]
            node_index = interested_nodes.index(node_label) + 1
            BendingXY_UR[node_index]['UR'] = [ur_1, ur_2, ur_3]
    odb3.close()      
    #endregion

    #region extract mass
    prop=mdb.models['Model-2'].rootAssembly.getMassProperties()
    mass['mass_wall_%s'%(j)]=prop['mass']/(length*width)
    mass['mass_face_%s'%(j)]=grammage_face[j]*10**(-12)
    mass['mass_core_%s'%(j)]=mass['mass_wall_%s'%(j)]-2*mass['mass_face_%s'%(j)]
    #endregion


    #region calculation
    #=============calculation===================

    #get the ABDH matrix of core only
    #A
    file_path = 'Model-1_elastic_properties(easycopy).txt'
    with open(file_path, 'r') as file:
        target_line_numbers = range(12)
        numbers = []
        for line_number, line in enumerate(file):
            if line_number in target_line_numbers:
                try:
                    number = float(line.strip())
                    numbers.append(number)
                except ValueError:
                    if line.strip() == "N/A":
                        numbers.append(0) 
                    else:
                        pass  
    file.close()

    Matrix['Q_matrix_wall_%s'%(j)]=np.array([ 
        [numbers[0]/(1-numbers[6]*numbers[8]), numbers[6]*numbers[1]/(1-numbers[6]*numbers[8]), 0],
        [numbers[6]*numbers[1]/(1-numbers[6]*numbers[8]), numbers[1]/(1-numbers[6]*numbers[8]), 0],
        [0, 0, numbers[3] ] ])

    Matrix['A_matrix_wall_%s'%(j)]=Matrix['Q_matrix_wall_%s'%(j)]*h[j]
    
    Matrix['Q_matrix_face_%s'%(j)]=np.array([ 
        [Mat_face[j][0]**2/(Mat_face[j][0]-Mat_face[j][2]**2*Mat_face[j][1]), Mat_face[j][0]*Mat_face[j][2]*Mat_face[j][1]/(Mat_face[j][0]-Mat_face[j][2]**2*Mat_face[j][1]), 0],
        [Mat_face[j][0]*Mat_face[j][2]*Mat_face[j][1]/(Mat_face[j][0]-Mat_face[j][2]**2*Mat_face[j][1]), Mat_face[j][0]*Mat_face[j][1]/(Mat_face[j][0]-Mat_face[j][2]**2*Mat_face[j][1]), 0],
        [0, 0, Mat_face[j][3]] ])
    
    Matrix['A_matrix_face_%s'%(j)]=Matrix['Q_matrix_face_%s'%(j)]*tf[j]
    
    Matrix['A_matrix_core_%s'%(j)]=Matrix['A_matrix_wall_%s'%(j)]-2*Matrix['A_matrix_face_%s'%(j)]
    #B=0
    #D
    jobname=['BendingX','BendingY','BendingXY']
    for i in jobname:
        #
        ur1u=globals()[i+'_UR'][1]['UR'][1]-globals()[i+'_UR'][5]['UR'][1]
        ur1d=globals()[i+'_UR'][3]['UR'][1]-globals()[i+'_UR'][7]['UR'][1]
        ur1=-(ur1u+ur1d)/2
        globals()['k11_'+i]=ur1*2/length
        #
        ur2u=globals()[i+'_UR'][6]['UR'][0]-globals()[i+'_UR'][5]['UR'][0]
        ur2d=globals()[i+'_UR'][8]['UR'][0]-globals()[i+'_UR'][7]['UR'][0]
        ur2=(ur2u+ur2d)/2
        globals()['k22_'+i]=ur2*2/width

        #
        ur12u=globals()[i+'_UR'][1]['UR'][0]-globals()[i+'_UR'][5]['UR'][0]
        ur12d=globals()[i+'_UR'][3]['UR'][0]-globals()[i+'_UR'][7]['UR'][0]
        ur12=(ur12u+ur12d)/2
        ur21u=globals()[i+'_UR'][6]['UR'][1]-globals()[i+'_UR'][5]['UR'][1]
        ur21d=globals()[i+'_UR'][8]['UR'][1]-globals()[i+'_UR'][7]['UR'][1]
        ur21=-(ur21u+ur21d)/2
        globals()['k12_'+i]=ur12*2/length+ur21*2/length

    import numpy as np
    
    curvature = np.array([
        [k11_BendingX, (k11_BendingY+k22_BendingX)/2, (k11_BendingXY+k12_BendingX)/2],
        [(k11_BendingY+k22_BendingX)/2, k22_BendingY, (k12_BendingY+k22_BendingXY)/2],
        [(k11_BendingXY+k12_BendingX)/2, (k12_BendingY+k22_BendingXY)/2, 2*k12_BendingXY]])
    
    Matrix['D_matrix_wall_%s'%(j)] = np.linalg.inv(curvature)

    Matrix['D_matrix_face_%s'%(j)]=Matrix['Q_matrix_face_%s'%(j)]*(tf[j]**3/12)

    Matrix['D_matrix_core_%s'%(j)]=Matrix['D_matrix_wall_%s'%(j)]-2*(Matrix['D_matrix_face_%s'%(j)]+(h[j]/2)**2*Matrix['A_matrix_face_%s'%(j)])

    #H
    Matrix['Qstar_matrix_wall_%s'%(j)]=np.array([ 
        [numbers[4], 0],
        [0,numbers[5]] ])

    Matrix['H_matrix_wall_%s'%(j)]=Matrix['Qstar_matrix_wall_%s'%(j)]*h[j]

    Matrix['Qstar_matrix_face_%s'%(j)]=np.array([ 
        [Mat_face[j][4], 0],
        [0,Mat_face[j][5]] ])

    Matrix['H_matrix_face_%s'%(j)]=Matrix['Qstar_matrix_face_%s'%(j)]*tf[j]

    Matrix['H_matrix_core_%s'%(j)]=Matrix['H_matrix_wall_%s'%(j)]-2*Matrix['H_matrix_face_%s'%(j)]

    mdb.models.changeKey(fromName='Model-1', toName='Tension-%s'%(j))
    mdb.models.changeKey(fromName='Model-2', toName='Bending-%s'%(j))
    mdb.Model(name='Model-1', modelType=STANDARD_EXPLICIT);
    mdb.Model(name='Model-2', modelType=STANDARD_EXPLICIT);

    #endregion

#region Total ABDH Matrix
# A&H of upper face
Matrix['Q_matrix_face_%s'%(wallnum)]=np.array([ 
    [Mat_face[wallnum][0]**2/(Mat_face[wallnum][0]-Mat_face[wallnum][2]**2*Mat_face[wallnum][1]), Mat_face[wallnum][0]*Mat_face[wallnum][2]*Mat_face[wallnum][1]/(Mat_face[wallnum][0]-Mat_face[wallnum][2]**2*Mat_face[wallnum][1]), 0],
    [Mat_face[wallnum][0]*Mat_face[wallnum][2]*Mat_face[wallnum][1]/(Mat_face[wallnum][0]-Mat_face[wallnum][2]**2*Mat_face[wallnum][1]), Mat_face[wallnum][0]*Mat_face[wallnum][1]/(Mat_face[wallnum][0]-Mat_face[wallnum][2]**2*Mat_face[wallnum][1]), 0],
    [0, 0, Mat_face[wallnum][3]] ])

A_matrix=Matrix['A_matrix_face_%s'%(wallnum)]=Matrix['Q_matrix_face_%s'%(wallnum)]*tf[wallnum]#initialization

Matrix['D_matrix_face_%s'%(wallnum)]=Matrix['Q_matrix_face_%s'%(wallnum)]*(tf[wallnum]**3/12)

Matrix['Qstar_matrix_face_%s'%(wallnum)]=np.array([ 
[Mat_face[wallnum][4], 0],
[0,Mat_face[wallnum][5]] ])

H_matrix=Matrix['H_matrix_face_%s'%(wallnum)]=Matrix['Qstar_matrix_face_%s'%(wallnum)]*tf[wallnum]#initialization

#=============AH-Matrix=====================
for k in range(wallnum):
    A_matrix=A_matrix+Matrix['A_matrix_face_%s'%(k)]+Matrix['A_matrix_core_%s'%(k)]
    H_matrix=H_matrix+Matrix['H_matrix_face_%s'%(k)]+Matrix['H_matrix_core_%s'%(k)]

#correction for H
if wallnum>=2:
    H_matrix[0][0]=H_matrix[0][0]*cH11**(wallnum-1)
    H_matrix[1][1]=H_matrix[1][1]*cH22**(wallnum-1)


#=============B-Matrix=====================
if wallnum ==1:
    B_matrix =  Matrix['A_matrix_face_%s'%(0)]*-(sum(h[:1])/2-0)+\
                Matrix['A_matrix_face_%s'%(1)]*-(sum(h[:1])/2-h[0])+\
                Matrix['A_matrix_core_%s'%(0)]*-(sum(h[:1])/2-h[0]/2)

elif wallnum ==2:
    B_matrix=   Matrix['A_matrix_face_%s'%(0)]*-(sum(h[:2])/2-0)+\
                Matrix['A_matrix_face_%s'%(1)]*-(sum(h[:2])/2-h[0])+\
                Matrix['A_matrix_face_%s'%(2)]*-(sum(h[:2])/2-h[1]-h[0])+\
                Matrix['A_matrix_core_%s'%(0)]*-(sum(h[:2])/2-h[0]/2)+\
                Matrix['A_matrix_core_%s'%(1)]*-(sum(h[:2])/2-h[0]-h[1]/2)
    
elif wallnum==3:
    B_matrix=   Matrix['A_matrix_face_%s'%(0)]*-(sum(h)/2-0)+\
                Matrix['A_matrix_face_%s'%(1)]*-(sum(h)/2-h[0])+\
                Matrix['A_matrix_face_%s'%(2)]*-(sum(h)/2-h[1]-h[0])+\
                Matrix['A_matrix_face_%s'%(3)]*-(sum(h)/2-h[2]-h[1]-h[0])+\
                Matrix['A_matrix_core_%s'%(0)]*-(sum(h)/2-h[0]/2)+\
                Matrix['A_matrix_core_%s'%(1)]*-(sum(h)/2-h[0]-h[1]/2)+\
                Matrix['A_matrix_core_%s'%(2)]*-(sum(h)/2-h[0]-h[1]-h[2]/2)    

#=============D-Matrix=====================
if wallnum ==1:
    D_matrix=\
            Matrix['D_matrix_face_%s'%(0)]+Matrix['A_matrix_face_%s'%(0)]*(0-sum(h[:1])/2)**2+\
            Matrix['D_matrix_face_%s'%(1)]+Matrix['A_matrix_face_%s'%(1)]*(sum(h[:1])-sum(h[:1])/2)**2+\
            Matrix['D_matrix_core_%s'%(0)]

elif wallnum ==2:
    D_matrix=\
            Matrix['D_matrix_face_%s'%(0)]+Matrix['A_matrix_face_%s'%(0)]*(0-sum(h[:2])/2)**2+\
            Matrix['D_matrix_face_%s'%(1)]+Matrix['A_matrix_face_%s'%(1)]*(sum(h[:1])-sum(h[:2])/2)**2+\
            Matrix['D_matrix_face_%s'%(2)]+Matrix['A_matrix_face_%s'%(2)]*(sum(h[:2])-sum(h[:2])/2)**2+\
            Matrix['D_matrix_core_%s'%(0)]+Matrix['A_matrix_core_%s'%(0)]*(h[0]/2-sum(h[:2])/2)**2+\
            Matrix['D_matrix_core_%s'%(1)]+Matrix['A_matrix_core_%s'%(1)]*(h[0]+h[1]/2-sum(h[:2])/2)**2
    
elif wallnum==3:
    D_matrix=\
            Matrix['D_matrix_face_%s'%(0)]+Matrix['A_matrix_face_%s'%(0)]*(0-sum(h)/2)**2+\
            Matrix['D_matrix_face_%s'%(1)]+Matrix['A_matrix_face_%s'%(1)]*(sum(h[:1])-sum(h)/2)**2+\
            Matrix['D_matrix_face_%s'%(2)]+Matrix['A_matrix_face_%s'%(2)]*(sum(h[:2])-sum(h)/2)**2+\
            Matrix['D_matrix_face_%s'%(3)]+Matrix['A_matrix_face_%s'%(3)]*(sum(h)-sum(h)/2)**2+\
            Matrix['D_matrix_core_%s'%(0)]+Matrix['A_matrix_core_%s'%(0)]*(h[0]/2-sum(h)/2)**2+\
            Matrix['D_matrix_core_%s'%(1)]+Matrix['A_matrix_core_%s'%(1)]*(h[0]+h[1]/2-sum(h)/2)**2+\
            Matrix['D_matrix_core_%s'%(2)]+Matrix['A_matrix_core_%s'%(2)]*(h[0]+h[1]+h[2]/2-sum(h)/2)**2



#check positive defined
top = np.hstack((A_matrix, B_matrix))
bottom = np.hstack((B_matrix, D_matrix))
ABD_matrix = np.vstack((top, bottom))

if not (np.allclose(ABD_matrix, ABD_matrix.T) and np.all(np.linalg.eigvals(ABD_matrix) > 0)):
    print("Error! ABDH_Matrix is not positive defined. Please check the code.")
    error()
if not (np.allclose(H_matrix, H_matrix.T) and np.all(np.linalg.eigvals(H_matrix) > 0)):
    print("Error! ABDH_Matrix is not positive defined. Please check the code.")
    error()


#endregion

#region total mass
density=mass['mass_face_%s'%(wallnum)]=grammage_face[wallnum]*10**(-12)
for k in range(wallnum):
    density=density+mass['mass_face_%s'%(k)]+mass['mass_core_%s'%(k)]

#endregion

#region output

#========print
print("--------ABD-Matrix Calculation Complete--------\n")
print("Units:mm-Mpa-t-N-s")
print("A=\n")
rows, cols = A_matrix.shape
for i in range(rows):
    row_str = " ".join("{:10.5f}".format(cell) for cell in A_matrix[i])
    print(row_str)
print("B=\n")
rows, cols = B_matrix.shape
for i in range(rows):
    row_str = " ".join("{:10.5f}".format(cell) for cell in B_matrix[i])
    print(row_str)
print("D=\n")
rows, cols = D_matrix.shape
for i in range(rows):
    row_str = " ".join("{:10.5f}".format(cell) for cell in D_matrix[i])
    print(row_str)
print("H=\n")
rows, cols = H_matrix.shape
for i in range(rows):
    row_str = " ".join("{:10.5f}".format(cell) for cell in H_matrix[i])
    print(row_str)



#===========Write file============
if wallnum==1:
    corType='Sigle-wall'
elif wallnum==2:
    corType='Double-wall'
elif wallnum==3:
    corType='Triple-wall'

with open("ABDH_Matrix.txt", "w") as f:
    f.write("=====================================================\n")
    f.write("Corrugated Cardboard to Shell %s\n"%(Version))
    f.write("=====================================================\n")
    f.write("Parameter Information:\n")
    f.write("Units:mm-Mpa-t-N-s\n")    
    f.write("Type:%s\n"%(corType))
    f.write("Material:\n")
    for i in range(2*wallnum+1):
        if (i+1)%2!=0:
            f.write("Liner%s:\n"%(i/2+1))
            f.write("\tgrammage:%s\t%s\n"%(grammage_face[i/2],'g/m^2'))
            f.write("\tt:%s\n"%(tf[i/2]))
            f.write("\tElastic constant:")
            for j in Mat_face[i/2]:
                f.write('\t{:6.1f}'.format(j))
            f.write("\n")
            continue
        elif (i+1)%2==0:
            f.write("Flute%s:\n"%((i-1)/2+1))
            f.write("\tgrammage:%s\t%s\n"%(grammage_core[(i-1)/2],'g/m^2'))
            f.write("\tt:%s\n"%(tc[(i-1)/2]))
            f.write("\tp:%s\n"%(period[(i-1)/2]))
            f.write("\th:%s\n"%(h[(i-1)/2]))    
            f.write("\tElastic constant:")
            for j in Mat_core[(i-1)/2]:
                f.write('\t{:6.1f}'.format(j))
            f.write("\n")
            continue

    f.write("A=\n")
    rows, cols = len(A_matrix), len(A_matrix[0])
    for i in range(rows):
        row_str = " ".join("{:10.5f}".format(cell) for cell in A_matrix[i])
        f.write(row_str + "\n")

    f.write("B=\n")
    rows, cols = len(B_matrix), len(B_matrix[0])
    for i in range(rows):
        row_str = " ".join("{:10.5f}".format(cell) for cell in B_matrix[i])
        f.write(row_str + "\n")
    
    f.write("D=\n")   
    rows, cols = len(D_matrix), len(D_matrix[0])
    for i in range(rows):
        row_str = " ".join("{:10.5f}".format(cell) for cell in D_matrix[i])
        f.write(row_str + "\n")
    
    f.write("H=\n") 
   
    rows, cols = len(H_matrix), len(H_matrix[0])
    for i in range(rows):
        row_str = " ".join("{:10.5f}".format(cell) for cell in H_matrix[i])
        f.write(row_str + "\n")
    if wallnum>=2:
        f.write("Correction coefficient: cH11=%s\tcH22=%s\n"%(cH11,cH22)) 
    f.write("Density(mass per unit surface area)=%s\n"%(density)) 
    f.write("!!Density is only a reference value due to overlapping of shell element!!") 

print("--------Matrices are written to files--------\n")
#endregion

#*****************end*******************
end_time = time.time()
execution_time = end_time - start_time
print("========== END ==========\n")
print("Execution time: {:.2f} seconds".format(execution_time))





















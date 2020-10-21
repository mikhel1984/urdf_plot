# Use matplotlib to show URDF robot structure
# 2020, Stanislav Mikhel

# parse
from xml.dom import minidom
from math import sin, cos
# plot
import matplotlib.pyplot as plt 
import mpl_toolkits.mplot3d as plt3d
import numpy as np

Ln = 5E-2    # joint axis size
Tl = 3E-2    # frame size
Clrs = {"revolute":"purple", "prismatic":"orange"}  # joint type colors
Frms = ["red","green","blue"]                       # X,Y,Z colors
  
def Tx(q):
  return np.array([[1,0,0,q],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
    
def Ty(q):
  return np.array([[1,0,0,0],[0,1,0,q],[0,0,1,0],[0,0,0,1]])
    
def Tz(q):
  return np.array([[1,0,0,0],[0,1,0,0],[0,0,1,q],[0,0,0,1]])
    
def Txyz(x,y,z): return np.array([[1,0,0,x],[0,1,0,y],[0,0,1,z],[0,0,0,1]])
    
def Rx(q):
  s = sin(q); c = cos(q)
  return np.array([[1,0,0,0],[0,c,-s,0],[0,s,c,0],[0,0,0,1]])
  
def Ry(q):
  s = sin(q); c = cos(q)  
  return np.array([[c,0,s,0],[0,1,0,0],[-s,0,c,0],[0,0,0,1]])
  
def Rz(q):
  s = sin(q); c = cos(q)
  return np.array([[c,-s,0,0],[s,c,0,0],[0,0,1,0],[0,0,0,1]])

H = [Rx, Ry, Rz, Tx, Ty, Tz] 
  
def Hrot(lst):
  r,p,w = lst
  return np.array([
    [cos(r)*cos(p), -sin(r)*cos(w)+cos(r)*sin(p)*sin(w), sin(r)*sin(w)+cos(r)*sin(p)*cos(w), 0],
    [sin(r)*cos(p), cos(r)*cos(w)+sin(r)*sin(p)*sin(w), -cos(r)*sin(w)+sin(r)*sin(p)*cos(w), 0],
    [-sin(p), cos(p)*sin(w), cos(p)*cos(w), 0],
    [0,0,0,1]]) 

# file processor
class Robot:

  def __init__(self,fname):
    '''Build robot description from URDF
       fname - URDF file'''
    # read urdf 
    dom = minidom.parse(fname)
    dom.normalize()
    self._getJointDict(dom)
    self._getLinkDict(dom)
    # processing
    self._sortJoints()
    self.fname = fname
    
  def _getJointDict(self,dom):
    '''Read joint info from XML file
       dom - XML object '''
    # extract information
    joints = dom.getElementsByTagName("joint")
    jointDict = {}
    for J in joints:
      name = J.getAttribute("name")
      d = {"fn":None,"axis":None}    # save information into dictionary 
      d["type"] = J.getAttribute("type")
      parent = J.getElementsByTagName("parent")
      d["parent"] = parent[0].getAttribute("link") if parent else None 
      child = J.getElementsByTagName("child")
      d["child"] = child[0].getAttribute("link") if child else None 
      origin = J.getElementsByTagName("origin")
      xyz = self._getList(origin[0].getAttribute("xyz")) if origin else None
      rpy = self._getList(origin[0].getAttribute("rpy")) if origin else None
      d["H"] = self._toMatrix(xyz,rpy)
      axis = J.getElementsByTagName("axis")
      if axis:
        i = axis[0].getAttribute("xyz").split().index('1')
        d["axis"] = i 
        if d["type"] == 'revolute':
          d["fn"] = (H[i],d["H"]) 
        elif d["type"] == 'prismatic':
          d["fn"] = (H[i+3],d["H"]) 
      jointDict[name] = d
    self.joints = jointDict

  def _getLinkDict(self,dom):
    '''Read link info from XML file
       dom - XML object '''
    links = dom.getElementsByTagName("link") 
    linkDict = {} 
    for L in links:
      name = L.getAttribute("name") 
      linkDict[name] = {"child": [], "parent": None, "pos":np.identity(4)}
    self.links = linkDict 

  def _toMatrix(self,xyz,rpy):
    '''(Private) Convert positions to matrix 
       xyz - Cartesian deflection
       rpy - rotation
       Return homogenous matrix'''
    m = np.identity(4)
    if xyz: m = m.dot(Txyz(*xyz))
    if rpy: m = m.dot(Hrot(rpy))
    return m
    
  def _getList(self,s):
    '''(Private) Convert string into list
       s - string 
       Return list of numbers''' 
    return [float(v) for v in s.split()] if s else None 

  def _sortJoints(self):
    '''Find parents for all links'''
    last = None
    # connect tree elements
    for nm in self.joints:
      jnt = self.joints[nm]
      lnk = self.links[jnt.get("parent")]
      lnk["child"].append(nm)
      lnk = self.links[jnt.get("child")]
      lnk["parent"] = nm
      last = nm 
    # find root
    while True:
      jnt = self.joints[last] 
      lnk = self.links[jnt.get("parent")] 
      last = lnk.get("parent")
      if last is None:
        self.root = jnt.get("parent")
        break

  def _update(self,qdct):
    '''(Private) Calculate joint transformation
       qdct - dictionary of joint angles'''
    for key in qdct:
      q = qdct[key]
      jnt = self.joints[key] 
      fns = jnt["fn"] 
      jnt["H"] = fns[1].dot(fns[0](q)) 

  def _getSegments(self,qdct):
    '''(Private) Find position for each joint
       qdct - dictionary of joint angles
       Return lists of nodes, lines, axes'''
    self._update(qdct)
    stack = [self.root]
    lines, nodes, axs, frame = [], [], [], []  # collect parts
    while stack: 
      # get link
      link = self.links[stack.pop()] 
      # adjusent joint position
      pos = link["pos"]
      p0 = pos[:3,3]
      jnm = link.get("parent")
      # get type and axis
      tp = None
      if jnm:
        tp = self.joints[jnm]["type"]
        if tp in Clrs:
          k = self.joints[jnm]["axis"]
          p1 = pos[:3,0]
          p2 = pos[:3,1]
          p3 = pos[:3,2] 
          Z = pos[:3,k]
          aa = (p0-Ln*Z, p0+Ln*Z)
          axs.append((aa,tp))
          frame.append( ((p0,p0+Tl*p1),(p0,p0+Tl*p2),(p0,p0+Tl*p3)) )
      nodes.append((p0,tp))
      # find edges
      for j in link["child"]:
        jnt = self.joints[j]
        child = jnt.get("child")
        if child is not None:
          # get position
          link = self.links[child] 
          link["pos"] = pos.dot(jnt["H"])  
          p1 = link["pos"][:3,3] 
          # save
          lines.append((p0,p1))
          stack.append(child)
    return lines, nodes, axs, frame
        
  def plot(self,qdct={}):
    '''Visualize robot structure
       qdct - dictionary of joint angles'''
    lines, nodes, axes, frame = self._getSegments(qdct) 
    # make figure 
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    plt.title(self.fname)
    # add nodes
    for n in nodes:
      nd,nt = n
      if nt == 'revolute':
        ax.scatter(nd[0],nd[1],nd[2], color=Clrs[nt], marker='o', s=3) 
      elif nt == 'prismatic':
        ax.scatter(nd[0],nd[1],nd[2], color=Clrs[nt], marker='s', s=3) 
      else:
        ax.scatter(nd[0],nd[1],nd[2], marker=',', s=0.3) 
    # add lines 
    for ln in lines:
      xs = (ln[0][0],ln[1][0])
      ys = (ln[0][1],ln[1][1])
      zs = (ln[0][2],ln[1][2]) 
      l = plt3d.art3d.Line3D(xs,ys,zs,color='black')
      ax.add_line(l) 
    # add axis
    for aa,at in axes:
      xs = (aa[0][0],aa[1][0])
      ys = (aa[0][1],aa[1][1])
      zs = (aa[0][2],aa[1][2])
      l = plt3d.art3d.Line3D(xs,ys,zs,color=Clrs[at],linestyle='--')
      ax.add_line(l)
    # add frames
    for ll in frame:
      for i,l in enumerate(ll):
        xs = (l[0][0],l[1][0])
        ys = (l[0][1],l[1][1])
        zs = (l[0][2],l[1][2])
        t = plt3d.art3d.Line3D(xs,ys,zs,color=Frms[i],linewidth=2)
        ax.add_line(t)
    plt.show()

# Example:
# 
# rbt = Robot("my_robot.urdf")
# qs = {"joint1":0.7, "joint3":-0.3}
# rbt.plot(qs)



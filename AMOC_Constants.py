import numpy as np

#earth radius, in meter
rds = 6.371e6

#ocean fraction
fraco = np.array([0.16267, 0.22222, 0.21069])

#atmospheric box latitudes, in radians (two latitudes for northern box, two for southern)
latia = np.array([-90.0, -30.0, 45.0, 90.0])*(np.pi/180.0)
#oceanic box boundary
latio = np.array([-60.0, -30.0, 45.0, 80.0])*(np.pi/180.0)

#sine of latitudes
sinlatia = np.sin(latia)
sinlata = (sinlatia[1:] + sinlatia[:-1])/2

#box mass centers, in radians
lata = np.arcsin(sinlata)

#area of each atmospheric box, in square meter
areaa = (2 * np.pi * rds**2)*np.diff(sinlatia)

#meridional distance between boxes, in meter
ydis = rds*np.diff(lata)

#perimeter of boundaries, in meter
perim = (2 * np.pi* rds)*np.cos(latia[1:-1])

#atmospheric heat capacities, in J/K
ca = areaa * (5300 * 1004.0)

#area of each oceanic box
areao = np.diff(np.sin(latio))*(80.0*np.pi/180.0*rds**2)
areao = np.array([*areao, areao[1]])


#height of each oceanic box, near-surface current and deep current
z1, z2 = 600, 4000
heighto = np.array([z2, z1, z2, z2-z1])

#volume of each oceanic box
vo = areao*heighto

#heat capacity of unit mass sea water, in J/(kg K)
cswt = 1025 * 4200

#oceanic heat capacities, in J/K
co = vo * cswt

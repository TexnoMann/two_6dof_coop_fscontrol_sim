from spatialmath import *
from roboticstoolbox.robot import *
import math

l1 = [0.163,    0.0]
l2 = [0.007,    0.425]
l3 = [0.0,      0.392]
l4 = [0.127,    0.0]
l5 = [0.1,      0.0]
l6 = [0.1,      0.0]

L = [
    RevoluteDH(d=l1[0], a=l1[1], alpha=-math.pi/2),
    RevoluteDH(d=l2[0], a=l2[1], alpha=         0),
    RevoluteDH(d=l3[0], a=l3[1], alpha=         0),
    RevoluteDH(d=l4[0], a=l4[1], alpha=-math.pi/2),
    RevoluteDH(d=l5[0], a=l5[1], alpha= math.pi/2),
    RevoluteDH(d=l6[0], a=l6[1], alpha=         0)
]
robot = DHRobot(L, name="ur5e")

class ur5e(DHRobot):
    def __init__(self):
        super().__init__(L, name="ur5e")


q = [0, -math.pi/3, math.pi/3 + math.pi/3, -math.pi/3 - math.pi/2, -math.pi/2, 0]
# q = [0.0, -0.8004778081347315, 1.6989733070614719, 3.8164067555811325, -1.5783361491636159, 0.0]
print(robot.fkine(q))
plt = robot.teach(q)
plt.hold()
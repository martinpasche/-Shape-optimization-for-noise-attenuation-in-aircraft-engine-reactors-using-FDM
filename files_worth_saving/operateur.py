global dx, dy
dx = 0.001
dy = 0.001

def d_x(i,j):
    return [[1/(2*dx),[i+1,j]],[-1/(2*dx),[i-1,j]]]
def d_y(i,j):
    return [[1/(2*dy),[i,j+1]],[-1/(2*dy),[i,j-1]]]

def d_2_x(i,j):
    return [[1/(dx*dx),[i+1,j]],[-2/(dx*dx),[i,j]],[1/(dx*dx),[i-1,j]]]
def d_2_y(i,j):
    return [[1/(dy*dy),[i,j+1]],[-2/(dy*dy),[i,j]],[1/(dy*dy),[i,j-1]]]
def lap(i,j):
        return [[1/(dx*dx),[i+1,j]],[-2/(dx*dx),[i,j]],[1/(dx*dx),[i-1,j]],\
    [1/(dy*dy),[i,j+1]],[-2/(dy*dy),[i,j]],[1/(dy*dy),[i,j-1]]]


# the code here is meant to be copy/pasted into a STAN program after it has been output
from sympy import symbols, expand, integrate, diff, solve
import json
import regex as re

OUTFILE='complicated_bspline_code.stan'

# a00,a10,a20,a30,a01,a11,a21,a31,a02,a12,a22,a32,a03,a13,a23,a33= symbols('a00 a10 a20 a30 a01 a11 a21 a31 a02 a12 a22 a32 a03 a13 a23 a33')

x,y = symbols('x y')
p = 0

for i in range(4):
    for j in range(4):
        p = p + y**i *x**j *symbols('a%s%s' % (i,j))

print('p')
print(p)

dx = diff(p,x)
dy = diff(p,y)
dxy = diff(dx,y)

dxx = diff(dx,x)
dyy = diff(dy,y)
dyxx = diff(dxx,y)
dxyy = diff(dyy,x)

dyyy = diff(dyy,y)
dxxx = diff(dxx,x)

# I admit that this is a bit hand-wavy, but for
# d2s I rotated z = x * y by 45 degrees to z = (x^2-y^2)/2
# the magnitudes of the 2nd derivatives are 1 and 1 at all points
# the mixed partial derivative is also 1 at each point
d2s = dxx**2 + dyy**2 + 2 * (dxy)**2

# for the 3rd derivative, I used z = x^2 * y (+ y^2 * x) with the same rotation
# the magnitudes of the 3rd derivatives are 3/sqrt(2) and 3/sqrt(2) compared to 2 for one of the partials
# without the second term y^2*x
# with the second term included the magnitude becomes 

d3s = dyyy**2 + dxxx**2 + 3*(dxyy)**2 + 3*(dyxx) ** 2

print('dx')
print(dx)
print('dy')
print(dy)
print('dxy')
print(dxy)

# evaluate

# f_YX fdy_YX fdx_YX fxy_YX
f_symbols = symbols([
    'f%s_%s%s' % (prefix,i,j)
    for i in range(2)
    for j in range(2)
    for prefix in [
            '','dy','dx','dxy'
    ]
])

a_symbols = symbols(['a%s%s' % (i,j) for i in range(4) for j in range(4)])

expressions = []
for i in range(2):
    for j in range(2):
        v = {'x':j,'y':i}
        # constant
        expressions.append(
            symbols('f_%s%s' % (i,j)) - p.subs(v)
        )
        # derivative y
        expressions.append(
            symbols('fdy_%s%s' % (i,j)) - dy.subs(v)
        )        

        # derivative x
        expressions.append(
            symbols('fdx_%s%s' % (i,j)) - dx.subs(v)
        )                

        # mixed derivative xy
        expressions.append(
            symbols('fdxy_%s%s' % (i,j)) - dxy.subs(v)
        )

print(expressions)
        
solution = solve(expressions,a_symbols,dict=True)[0]
print(solution)

# now to metaprogram

lines = []

# define the blocks

mname = 'map'
mname_dy = 'map_1d_y'
mname_dx = 'map_1d_x'
mname_dxy = 'map_2d_xy'

prefix_pairs = [
    (mname, ''),
    (mname_dy,'dy'),
    (mname_dx,'dx'),
    (mname_dxy,'dxy'),
]

# block-shifted definitions

for i in range(2):
    for j in range(2):
        for prefix_pair in prefix_pairs:
            lines.append(
                f'matrix[ydim-1,xdim-1] f{prefix_pair[1]}_{i}{j}'
            )
            
for i in range(2):
    for j in range(2):
        for prefix_pair in prefix_pairs:
            lines.append(
                f'f{prefix_pair[1]}_{i}{j} = block({prefix_pair[0]},{i+1},{j+1},ydim-1,xdim-1)',
            )

#aYX definitions

for i in range(4):
    for j in range(4):
        lines.append(
            f'matrix[ydim-1,xdim-1] a{i}{j}'
        )

for i in range(4):
    for j in range(4):
        v = f'a{i}{j}'
        s = solution[symbols(v)]
        lines.append(
            f'{v}={s}'
        )


# penalty

ss = expand(dx**2 + dy**2)
integral = ss.integrate((y,0,1)).integrate((x,0,1))
istring = str(integral)

print(istring)

asterisk_pattern = r'(?<=a\d\d)\*(?=a\d\d)'
asterisk_replace = r'.*'

asterisk_pattern2 = r'(a..)\^2'
asterisk_replace2 = r'\1.*\1'
lines.append(
    'penalty += W_derivative * sum(%s)' %
    re.sub(asterisk_pattern2, asterisk_replace2, re.sub(asterisk_pattern,asterisk_replace,istring.replace('**','^')))
)

def add_indices(x):
    return re.sub('(a..)',r'\1[i,j]',x)

# write out expanded forms of p, dx, dy, dxy
lines.extend([
    add_indices(str(expand(p)).replace('**','^')),
    add_indices(str(expand(dy)).replace('**','^')),
    add_indices(str(expand(dx)).replace('**','^')),
    add_indices(str(expand(dxy)).replace('**','^'))
])


# d2s
integral2 = d2s.integrate((y,0,1)).integrate((x,0,1))
istring2 = str(integral2)
lines.append(
    'penalty += W_derivative2 * sum(%s)' %
    re.sub(asterisk_pattern2, asterisk_replace2, re.sub(asterisk_pattern,asterisk_replace,istring2.replace('**','^')))
)

# d3s
integral3 = d3s.integrate((y,0,1)).integrate((x,0,1))
istring3 = str(integral3)
lines.append(
    'penalty += W_derivative3 * sum(%s)' %
    re.sub(asterisk_pattern2, asterisk_replace2, re.sub(asterisk_pattern,asterisk_replace,istring3.replace('**','^')))
)

print(json.dumps(lines, indent=2))

with open(OUTFILE,'w') as f:
    f.write('\n'.join([l+';' for l in lines]))

print('wrote lines')




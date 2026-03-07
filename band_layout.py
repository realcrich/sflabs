"""
band_layout.py â€” v7: Inner ring / courtyard + BOH/Circ rooms

Grid legend:
  S=studio  1=1BR  2=2BR  3=3BR  K=corner
  B=BOH  C=Circ  O=courtyard(open)
  .=corridor  #=core  ~=parcel_area  (space)=outside
  !=overlap  ?=gap(inside building, unclaimed)

Run: .venv/bin/python band_layout.py

$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$   CHANGE LOG    $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    03/07/26 - CR: added documentation to defined functions (max_inscribed_rect(), 
                   compute_layout(), and rasterize()) and improved general readability
    ...
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
"""
import json, math
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPoly

# â”€â”€ Geometry â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
def poly_area(pts):
    n=len(pts); return abs(sum(pts[i][0]*pts[(i+1)%n][1]-pts[(i+1)%n][0]*pts[i][1] for i in range(n)))/2

def centroid(pts):
    return (sum(p[0] for p in pts)/len(pts), sum(p[1] for p in pts)/len(pts))

def scale_poly(pts,f):
    cx,cy=centroid(pts); return [(cx+(x-cx)*f,cy+(y-cy)*f) for x,y in pts]

def scale_to_area(pts,target):
    return scale_poly(pts,math.sqrt(target/max(poly_area(pts),1)))

def center_at_origin(pts):
    cx,cy=centroid(pts); return [(x-cx,y-cy) for x,y in pts]

def dist(a,b): return math.hypot(b[0]-a[0],b[1]-a[1])
def vadd(a,b): return (a[0]+b[0],a[1]+b[1])
def vscale(v,s): return (v[0]*s,v[1]*s)
def vneg(v): return (-v[0],-v[1])

def point_in_poly(pt,poly):
    x,y=pt;n=len(poly);inside=False;j=n-1
    for i in range(n):
        xi,yi=poly[i];xj,yj=poly[j]
        if ((yi>y)!=(yj>y)) and (x<(xj-xi)*(y-yi)/(yj-yi)+xi): inside=not inside
        j=i
    return inside

def line_intersect(p1,p2,p3,p4):
    d=(p1[0]-p2[0])*(p3[1]-p4[1])-(p1[1]-p2[1])*(p3[0]-p4[0])
    if abs(d)<1e-10: return None
    t=((p1[0]-p3[0])*(p3[1]-p4[1])-(p1[1]-p3[1])*(p3[0]-p4[0]))/d
    return (p1[0]+t*(p2[0]-p1[0]),p1[1]+t*(p2[1]-p1[1]))

def signed_area(pts):
    n=len(pts); return sum(pts[i][0]*pts[(i+1)%n][1]-pts[(i+1)%n][0]*pts[i][1] for i in range(n))/2

def inset_polygon(poly,d):
    n=len(poly);ccw=signed_area(poly)>0;lines=[]
    for i in range(n):
        j=(i+1)%n;dx,dy=poly[j][0]-poly[i][0],poly[j][1]-poly[i][1];ln=math.hypot(dx,dy)
        if ln<1e-10: continue
        if ccw: nx,ny=-dy/ln,dx/ln
        else: nx,ny=dy/ln,-dx/ln
        lines.append(((poly[i][0]+nx*d,poly[i][1]+ny*d),(poly[j][0]+nx*d,poly[j][1]+ny*d)))
    result=[]
    for i in range(len(lines)):
        j=(i+1)%len(lines);pt=line_intersect(lines[i][0],lines[i][1],lines[j][0],lines[j][1])
        if pt: result.append(pt)
    return result if len(result)>=3 else scale_poly(poly,0.8)
'''
# â”€â”€ Max inscribed rectangle â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
def max_inscribed_rect(poly):
    best_area=0;best=None;angles=set();n=len(poly)
    for i in range(n):
        j=(i+1)%n;a=math.atan2(poly[j][1]-poly[i][1],poly[j][0]-poly[i][0])
        angles.add(a);angles.add(a+math.pi/2)
    for k in range(36): angles.add(math.pi*k/36)
    for angle in angles:
        ca,sa=math.cos(-angle),math.sin(-angle)
        rot=[(x*ca-y*sa,x*sa+y*ca) for x,y in poly]
        all_y=sorted(p[1] for p in rot);min_y,max_y=all_y[0],all_y[-1]
        ns=50;dy=(max_y-min_y)/ns;levels=[]
        for si in range(ns+1):
            y=min_y+dy*si;xs=[]
            for i in range(len(rot)):
                j=(i+1)%len(rot);y1,y2=rot[i][1],rot[j][1]
                if abs(y2-y1)>1e-10 and (y1-y)*(y2-y)<=0:
                    t=(y-y1)/(y2-y1);xs.append(rot[i][0]+t*(rot[j][0]-rot[i][0]))
            if len(xs)>=2: levels.append((y,min(xs),max(xs)))
        for i in range(len(levels)):
            xl,xr=levels[i][1],levels[i][2]
            for j in range(i+1,len(levels)):
                xl=max(xl,levels[j][1]);xr=min(xr,levels[j][2])
                w=xr-xl;h=levels[j][0]-levels[i][0]
                if w>0 and h>0 and w*h>best_area:
                    best_area=w*h;rcx,rcy=(xl+xr)/2,(levels[i][0]+levels[j][0])/2
                    ca2,sa2=math.cos(angle),math.sin(angle)
                    best=(rcx*ca2-rcy*sa2,rcx*sa2+rcy*ca2,w,h,angle)
    return best
'''

def max_inscribed_rect(poly):
    """
    Finds the maximum-area rectangle that fits inside a given polygon.

    The approach works by rotating the polygon to many candidate angles, then
    using horizontal scanlines to find the widest horizontal span at each level.
    It then tests all pairs of scanlines to find the largest axis-aligned
    rectangle in the rotated frame, before rotating the result back.

    Args:
        poly: A list of (x, y) tuples representing the polygon's vertices,
              in order (clockwise or counter-clockwise).

    Returns:
        A tuple of (cx, cy, width, height, angle) where:
            - (cx, cy) is the rectangle's center in the original coordinate space
            - width and height are the rectangle's dimensions
            - angle is the rotation angle (in radians) of the rectangle
        Returns None if no valid rectangle is found.
    """

    best_area = 0
    best = None

    # --- Step 1: Build a set of candidate rotation angles to try ---
    # The optimal rectangle is often aligned with one of the polygon's edges,
    # so we add each edge angle (and its perpendicular) as a candidate.
    angles = set()
    n = len(poly)

    for i in range(n):
        j = (i + 1) % n
        edge_angle = math.atan2(
            poly[j][1] - poly[i][1],
            poly[j][0] - poly[i][0]
        )
        angles.add(edge_angle)
        angles.add(edge_angle + math.pi / 2)  # Perpendicular to the edge

    # Also add 36 evenly-spaced angles (every 5 degrees) as a general fallback
    for k in range(36):
        angles.add(math.pi * k / 36)

    # --- Step 2: For each candidate angle, rotate and scan ---
    for angle in angles:

        # Rotate the polygon by -angle so the candidate rectangle
        # becomes axis-aligned, making it easier to measure
        ca, sa = math.cos(-angle), math.sin(-angle)
        rotated_poly = [(x * ca - y * sa, x * sa + y * ca) for x, y in poly]

        # Find the vertical extent of the rotated polygon
        all_y = sorted(p[1] for p in rotated_poly)
        min_y, max_y = all_y[0], all_y[-1]

        # --- Step 3: Cast horizontal scanlines through the rotated polygon ---
        # Divide the vertical range into `num_slices` evenly-spaced levels
        num_slices = 50
        dy = (max_y - min_y) / num_slices
        levels = []  # Each entry: (y, left_x, right_x) — the horizontal span at that y

        for si in range(num_slices + 1):
            y = min_y + dy * si
            x_intersections = []

            # Find where the horizontal line y intersects the polygon edges
            for i in range(len(rotated_poly)):
                j = (i + 1) % len(rotated_poly)
                y1, y2 = rotated_poly[i][1], rotated_poly[j][1]

                # Only consider edges that cross this y level (skip horizontal edges)
                edge_crosses_y = abs(y2 - y1) > 1e-10 and (y1 - y) * (y2 - y) <= 0
                if edge_crosses_y:
                    t = (y - y1) / (y2 - y1)  # Interpolation parameter
                    x_intersect = rotated_poly[i][0] + t * (rotated_poly[j][0] - rotated_poly[i][0])
                    x_intersections.append(x_intersect)

            # We need at least two intersections to define a horizontal span
            if len(x_intersections) >= 2:
                levels.append((y, min(x_intersections), max(x_intersections)))

        # --- Step 4: Find the best rectangle from pairs of scanlines ---
        # For each pair of scanlines (top and bottom of the rectangle),
        # the valid horizontal span is the intersection of all spans between them.
        # This ensures the rectangle stays fully inside the polygon.
        for i in range(len(levels)):
            xl, xr = levels[i][1], levels[i][2]  # Start with the top scanline's span

            for j in range(i + 1, len(levels)):
                # Tighten the horizontal span to stay within all intermediate scanlines
                xl = max(xl, levels[j][1])
                xr = min(xr, levels[j][2])

                w = xr - xl
                h = levels[j][0] - levels[i][0]

                if w > 0 and h > 0 and w * h > best_area:
                    best_area = w * h

                    # Compute the rectangle's center in the rotated frame
                    rcx = (xl + xr) / 2
                    rcy = (levels[i][0] + levels[j][0]) / 2

                    # Rotate the center back to the original coordinate space
                    ca2, sa2 = math.cos(angle), math.sin(angle)
                    cx = rcx * ca2 - rcy * sa2
                    cy = rcx * sa2 + rcy * ca2

                    best = (cx, cy, w, h, angle)

    return best

def get_rect_corners(cx,cy,w,h,angle):
    ca,sa=math.cos(angle),math.sin(angle);hw,hh=w/2,h/2
    return [(cx+(-hw)*ca-(-hh)*sa,cy+(-hw)*sa+(-hh)*ca),(cx+(hw)*ca-(-hh)*sa,cy+(hw)*sa+(-hh)*ca),
            (cx+(hw)*ca-(hh)*sa,cy+(hw)*sa+(hh)*ca),(cx+(-hw)*ca-(hh)*sa,cy+(-hw)*sa+(hh)*ca)]

# â”€â”€ Parcels â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
def geo2ft(coords):
    ring=coords[:-1] if (len(coords)>1 and coords[-1]==coords[0]) else coords
    clng=sum(c[0] for c in ring)/len(ring);clat=sum(c[1] for c in ring)/len(ring)
    return [((lng-clng)*288200,(lat-clat)*364000) for lng,lat in ring]

PARCELS = {
    'p1': geo2ft([[-118.37490349962798,34.165693990658085],[-118.37467334233435,34.165789208630216],[-118.37469682777262,34.1659368934356],[-118.37484948312057,34.166181738727914],[-118.37480486078779,34.16620700050339],[-118.37496925885486,34.166504311600875],[-118.3753004035324,34.16637023025507]]),
    'p4': geo2ft([[-118.47128304946774,34.026451245666564],[-118.47099221555277,34.0262102081888],[-118.47056275986552,34.0265706378093],[-118.47063886593685,34.02664497622831],[-118.47054101527391,34.0267508520451],[-118.47073671659979,34.02690403384527]]),
    'p7': geo2ft([[-118.59721278909241,34.1881917698952],[-118.59721614524025,34.18744774676111],[-118.59636368366544,34.18745052297892],[-118.5963603275176,34.18819732228182]]),
    'p9': geo2ft([[-118.25022093172518,34.0502674231122],[-118.24971135075477,34.049947450481866],[-118.24956967353441,34.05010270362472],[-118.24968164424081,34.050168970123906],[-118.24962908656231,34.050214409978864],[-118.2496610781927,34.0502352365705],[-118.24943942189637,34.05047379535168],[-118.24981875122825,34.05071424677597]]),
}

UNIT_CHAR = {'studio':'S','1br':'1','2br':'2','3br':'3','corner':'K'}
COLORS = {'studio':'#c8e6c9','1br':'#bbdefb','2br':'#ffe0b2','3br':'#f8bbd0',
          'corner':'#e0e0e0','boh':'#ffccbc','circ':'#d1c4e9'}
LABELS = {'studio':'ST','1br':'1B','2br':'2B','3br':'3B','corner':'CR'}

'''
# â”€â”€ Layout engine â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
def compute_layout(data):
    pid = data['project_id']
    fp_area = data['building']['floor_plate_sf']
    corr_w = data['circulation']['corridor_width_ft']
    d_floors = data['building']['story_distribution']['dwelling']
    avg_depth = sum(u['depth_ft'] for u in data['dwelling_units']) / len(data['dwelling_units'])

    parcel = PARCELS.get(pid)
    boundary = center_at_origin(scale_to_area(parcel, fp_area)) if parcel else \
               (lambda h:[(-h,-h),(h,-h),(h,h),(-h,h)])(math.sqrt(fp_area)/2)

    mir = max_inscribed_rect(boundary)
    if not mir:
        h=math.sqrt(fp_area)/2*0.9; mir=(0,0,h*2,h*2,0)
    rcx,rcy,rw,rh,rangle = mir
    rect = get_rect_corners(rcx,rcy,rw,rh,rangle)

    depth = avg_depth
    min_dim = min(rw,rh)
    max_depth = (min_dim - corr_w - 8) / 2
    if depth > max_depth: depth = max(12, max_depth)

    door_w = 4
    cs = depth + door_w

    ca,sa = math.cos(rangle),math.sin(rangle)
    dir_x,dir_y = (ca,sa),(-sa,ca)
    edge_dirs = [dir_x,dir_y,vneg(dir_x),vneg(dir_y)]
    edge_normals = [dir_y,vneg(dir_x),vneg(dir_y),dir_x]
    edge_lens = [rw,rh,rw,rh]

    # â”€â”€ Outer corner units â”€â”€
    corners = []
    for ci in range(4):
        v = rect[ci]
        o0 = v
        o1 = vadd(v, vscale(edge_dirs[ci], cs))
        i0 = vadd(o0, vscale(edge_normals[ci], depth))
        i1 = vadd(o1, vscale(edge_normals[ci], depth))
        corners.append(([o0, o1, i1, i0], 'corner'))

    # Unit queue â€” interleaved types
    queue=[]
    type_counts=[(u['type'],u['width_ft'],max(1,round(u['count']/d_floors))) for u in data['dwelling_units']]
    max_c=max(c for _,_,c in type_counts)
    for r in range(max_c):
        for utype,w,c in type_counts:
            if r<c: queue.append((utype,w))

    # â”€â”€ Outer regular units â”€â”€
    units=[];qi=0;units_by_edge={}
    for ei in range(4):
        edge_start=len(units)
        avail=edge_lens[ei]-cs
        if avail<10:
            units_by_edge[ei]=[]
            continue
        origin=vadd(rect[ei],vscale(edge_dirs[ei],cs))
        uw=queue[qi%len(queue)][1] if queue else 20
        num=max(1,round(avail/uw))
        actual_w=avail/num
        while actual_w<14 and num>1: num-=1; actual_w=avail/num
        for k in range(num):
            o0=vadd(origin,vscale(edge_dirs[ei],k*actual_w))
            o1=vadd(origin,vscale(edge_dirs[ei],(k+1)*actual_w))
            i0=vadd(o0,vscale(edge_normals[ei],depth))
            i1=vadd(o1,vscale(edge_normals[ei],depth))
            utype=queue[qi%len(queue)][0];qi+=1
            units.append(([o0,o1,i1,i0],utype))
        units_by_edge[ei]=list(range(edge_start,len(units)))

    # â”€â”€ Core geometry â”€â”€
    inner_w,inner_h=rw-2*depth,rh-2*depth
    corr_rect=get_rect_corners(rcx,rcy,inner_w,inner_h,rangle) if inner_w>0 and inner_h>0 else None
    core_w,core_h=inner_w-2*corr_w,inner_h-2*corr_w

    raw_core = get_rect_corners(rcx,rcy,core_w,core_h,rangle) if core_w>0 and core_h>0 else None

    # Chamfered core for rendering
    if core_w>0 and core_h>0:
        core_chamfer=min(corr_w*1.5,core_w/4,core_h/4)
        core_rect=[]
        for ci2 in range(4):
            v=raw_core[ci2]
            nxt=raw_core[(ci2+1)%4]; prv=raw_core[(ci2-1)%4]
            dn=dist(v,nxt); dp=dist(v,prv)
            dn_dir=((nxt[0]-v[0])/dn,(nxt[1]-v[1])/dn)
            dp_dir=((prv[0]-v[0])/dp,(prv[1]-v[1])/dp)
            core_rect.append(vadd(v,vscale(dp_dir,core_chamfer)))
            core_rect.append(vadd(v,vscale(dn_dir,core_chamfer)))
    else:
        core_rect=None

    # â”€â”€ Inner layout mode â”€â”€
    inner_units = []
    inner_corners = []
    courtyard_rect = None
    boh_room = None
    circ_room = None
    layout_mode = 'basic'
    inner_depth = 0

    min_core = min(core_w, core_h) if (core_w > 0 and core_h > 0) else 0
    if min_core > 0:
        # Inner depth may be smaller than outer to fit courtyard
        inner_depth = min(depth, (min_core - 2*corr_w - 8) / 2)

    if inner_depth >= 12 and core_w > 0 and core_h > 0:
        # â”€â”€ Case 1: Courtyard â€” inner ring of units â”€â”€
        layout_mode = 'courtyard'
        inner_cs = min(inner_depth + door_w, min(core_w, core_h) / 2)

        for ci in range(4):
            v = raw_core[ci]
            o0 = v
            o1 = vadd(v, vscale(edge_dirs[ci], inner_cs))
            i0 = vadd(o0, vscale(edge_normals[ci], inner_depth))
            i1 = vadd(o1, vscale(edge_normals[ci], inner_depth))
            inner_corners.append(([o0, o1, i1, i0], 'corner'))

        inner_edge_lens = [core_w, core_h, core_w, core_h]
        inner_units_by_edge = {}
        for ei in range(4):
            edge_start = len(inner_units)
            avail = inner_edge_lens[ei] - inner_cs
            if avail < 10:
                inner_units_by_edge[ei] = []
                continue
            origin = vadd(raw_core[ei], vscale(edge_dirs[ei], inner_cs))
            uw = queue[qi % len(queue)][1] if queue else 20
            num = max(1, round(avail / uw))
            actual_w = avail / num
            while actual_w < 14 and num > 1: num -= 1; actual_w = avail / num
            for k in range(num):
                o0 = vadd(origin, vscale(edge_dirs[ei], k * actual_w))
                o1 = vadd(origin, vscale(edge_dirs[ei], (k + 1) * actual_w))
                i0 = vadd(o0, vscale(edge_normals[ei], inner_depth))
                i1 = vadd(o1, vscale(edge_normals[ei], inner_depth))
                utype = queue[qi % len(queue)][0]; qi += 1
                inner_units.append(([o0, o1, i1, i0], utype))
            inner_units_by_edge[ei] = list(range(edge_start, len(inner_units)))

        cy_w = core_w - 2*inner_depth - 2*corr_w
        cy_h = core_h - 2*inner_depth - 2*corr_w
        courtyard_rect = get_rect_corners(rcx, rcy, cy_w, cy_h, rangle) if cy_w > 0 and cy_h > 0 else None

    elif core_w > 10 and core_h > 8:
        # â”€â”€ Case 2: BOH + Circ rooms in core â”€â”€
        layout_mode = 'core_rooms'
        boh_w, boh_h = 4, 4
        circ_w, circ_h = 6, 4
        total_w = boh_w + circ_w  # 10
        if core_w >= core_h:
            boh_dx = -total_w/2 + boh_w/2
            circ_dx = total_w/2 - circ_w/2
            boh_room = get_rect_corners(rcx + boh_dx*ca, rcy + boh_dx*sa, boh_w, boh_h, rangle)
            circ_room = get_rect_corners(rcx + circ_dx*ca, rcy + circ_dx*sa, circ_w, circ_h, rangle)
        else:
            boh_dy = -total_w/2 + boh_w/2
            circ_dy = total_w/2 - circ_w/2
            boh_room = get_rect_corners(rcx - boh_dy*sa, rcy + boh_dy*ca, boh_h, boh_w, rangle)
            circ_room = get_rect_corners(rcx - circ_dy*sa, rcy + circ_dy*ca, circ_h, circ_w, rangle)

    elif core_w > 4 and core_h > 4:
        # â”€â”€ Case 3: Notch â€” small rooms in core â”€â”€
        layout_mode = 'notch'
        boh_w, boh_h = 4, 4
        circ_w, circ_h = 6, 4
        total_w = boh_w + circ_w
        if core_w >= 10:
            boh_dx = -total_w/2 + boh_w/2
            circ_dx = total_w/2 - circ_w/2
            boh_room = get_rect_corners(rcx + boh_dx*ca, rcy + boh_dx*sa,
                                        boh_w, min(boh_h, core_h*0.9), rangle)
            circ_room = get_rect_corners(rcx + circ_dx*ca, rcy + circ_dx*sa,
                                         circ_w, min(circ_h, core_h*0.9), rangle)
        elif core_h >= 10:
            boh_dy = -total_w/2 + boh_w/2
            circ_dy = total_w/2 - circ_w/2
            boh_room = get_rect_corners(rcx - boh_dy*sa, rcy + boh_dy*ca,
                                        min(boh_h, core_w*0.9), boh_w, rangle)
            circ_room = get_rect_corners(rcx - circ_dy*sa, rcy + circ_dy*ca,
                                         min(circ_h, core_w*0.9), circ_w, rangle)
        else:
            boh_room = get_rect_corners(rcx, rcy, min(4, core_w*0.9), min(4, core_h*0.9), rangle)

    # â”€â”€ Always place BOH + Circ â”€â”€
    # For core_rooms they're in the core already.
    # For courtyard: on outer edge of inner ring units (corridor-facing), split on long edges.
    # For basic: on inner edge of outer ring units (corridor-facing), split on long edges.
    if boh_room is None:
        boh_w, boh_h = 4, 4
        circ_w, circ_h = 6, 4

        def _room_on_unit(unit_verts, rw, rh, on_outer):
            """Place room rect on the outer or inner edge of a unit, centered."""
            o0, o1, i1, i0 = unit_verts
            if on_outer:
                p0, p1, into = o0, o1, i0
            else:
                p0, p1, into = i0, i1, o0
            e_len = dist(p0, p1)
            e_dir = ((p1[0]-p0[0])/e_len, (p1[1]-p0[1])/e_len)
            d_len = dist(p0, into)
            d_dir = ((into[0]-p0[0])/d_len, (into[1]-p0[1])/d_len)
            off = max(0, (e_len - rw) / 2)
            s = vadd(p0, vscale(e_dir, off))
            return [s, vadd(s, vscale(e_dir, rw)),
                    vadd(vadd(s, vscale(e_dir, rw)), vscale(d_dir, rh)),
                    vadd(s, vscale(d_dir, rh))]

        if layout_mode == 'courtyard' and inner_units:
            # Outer edge of inner ring units, split on the two longest edges
            elens = [core_w, core_h, core_w, core_h]
            ranked = sorted(range(4), key=lambda e: len(inner_units_by_edge.get(e,[])), reverse=True)
            long_edges = [e for e in ranked if inner_units_by_edge.get(e)][:2]
            for idx, (rw2, rh2) in enumerate([(boh_w, boh_h), (circ_w, circ_h)]):
                ei = long_edges[idx % len(long_edges)]
                ulist = inner_units_by_edge[ei]
                ui = ulist[len(ulist) // 2]
                room = _room_on_unit(inner_units[ui][0], rw2, rh2, on_outer=True)
                if idx == 0: boh_room = room
                else: circ_room = room
        elif units:
            # Inner edge of outer ring units, split on the two longest edges
            ranked = sorted(range(4), key=lambda e: len(units_by_edge.get(e,[])), reverse=True)
            long_edges = [e for e in ranked if units_by_edge.get(e)][:2]
            for idx, (rw2, rh2) in enumerate([(boh_w, boh_h), (circ_w, circ_h)]):
                ei = long_edges[idx % len(long_edges)]
                ulist = units_by_edge[ei]
                ui = ulist[len(ulist) // 2]
                room = _room_on_unit(units[ui][0], rw2, rh2, on_outer=False)
                if idx == 0: boh_room = room
                else: circ_room = room

    return {
        'boundary':boundary,'rect':rect,'corr_rect':corr_rect,'core_rect':core_rect,
        'units':units,'corners':corners,'depth':depth,'corr_w':corr_w,
        'fp_area':fp_area,'pid':pid,'rw':rw,'rh':rh,'avg_depth':avg_depth,
        'layout_mode':layout_mode,
        'inner_units':inner_units,'inner_corners':inner_corners,
        'courtyard_rect':courtyard_rect,
        'boh_room':boh_room,'circ_room':circ_room,
        'core_w':core_w,'core_h':core_h,'inner_depth':inner_depth,
        'rcx':rcx,'rcy':rcy,'rangle':rangle,
    }
'''

def compute_layout(data):
    """
    Computes the full 2D floor plan layout for a residential building floor plate.

    Given a project's configuration data, this function:
      1. Retrieves or synthesizes the parcel boundary and fits a maximum inscribed rectangle.
      2. Places outer corner and regular dwelling units around the perimeter.
      3. Determines the core geometry (corridor, service core) inside the unit ring.
      4. Selects a layout mode for the core based on available space:
           - 'courtyard'   : inner ring of units around an open courtyard
           - 'core_rooms'  : back-of-house (BOH) and circulation rooms in the core
           - 'notch'       : smaller rooms squeezed into a tight core
           - 'basic'       : no core rooms, falls back to BOH/circ on unit edges
      5. Always ensures BOH and circulation rooms are placed somewhere.

    Args:
        data: A dict with keys:
            - 'project_id'        : str, used to look up parcel geometry in PARCELS
            - 'building'          : dict with 'floor_plate_sf' and
                                    'story_distribution' -> {'dwelling': int}
            - 'circulation'       : dict with 'corridor_width_ft'
            - 'dwelling_units'    : list of dicts with 'type', 'width_ft', 'depth_ft', 'count'

    Returns:
        A dict containing all computed geometry and metadata:
            boundary        - polygon of the floor plate boundary
            rect            - 4 corners of the max inscribed rectangle
            corr_rect       - 4 corners of the corridor rectangle (inside unit ring)
            core_rect       - chamfered polygon of the service core
            units           - list of (polygon, type) for outer regular units
            corners         - list of (polygon, 'corner') for outer corner units
            depth           - actual unit depth used (ft)
            corr_w          - corridor width (ft)
            fp_area         - floor plate area (sf)
            pid             - project ID
            rw, rh          - width and height of the inscribed rectangle
            avg_depth       - average unit depth from input data
            layout_mode     - one of 'basic', 'courtyard', 'core_rooms', 'notch'
            inner_units     - list of (polygon, type) for courtyard-ring units
            inner_corners   - list of (polygon, 'corner') for courtyard-ring corners
            courtyard_rect  - polygon of the open courtyard (or None)
            boh_room        - polygon of the back-of-house room (or None)
            circ_room       - polygon of the circulation room (or None)
            core_w, core_h  - dimensions of the service core
            inner_depth     - unit depth used for the inner ring
            rcx, rcy        - center of the inscribed rectangle
            rangle          - rotation angle of the inscribed rectangle (radians)
    """

    # -------------------------------------------------------------------------
    # Step 1: Extract inputs
    # -------------------------------------------------------------------------
    pid = data['project_id']
    fp_area = data['building']['floor_plate_sf']
    corr_w = data['circulation']['corridor_width_ft']
    d_floors = data['building']['story_distribution']['dwelling']
    avg_depth = sum(u['depth_ft'] for u in data['dwelling_units']) / len(data['dwelling_units'])

    # -------------------------------------------------------------------------
    # Step 2: Define the floor plate boundary
    # Use the parcel shape if available; otherwise synthesize a square.
    # -------------------------------------------------------------------------
    parcel = PARCELS.get(pid)
    if parcel:
        boundary = center_at_origin(scale_to_area(parcel, fp_area))
    else:
        h = math.sqrt(fp_area) / 2
        boundary = [(-h, -h), (h, -h), (h, h), (-h, h)]

    # -------------------------------------------------------------------------
    # Step 3: Find the maximum inscribed rectangle within the boundary
    # -------------------------------------------------------------------------
    mir = max_inscribed_rect(boundary)
    if not mir:
        # Fallback: use a slightly shrunk square centered at origin
        h = math.sqrt(fp_area) / 2 * 0.9
        mir = (0, 0, h * 2, h * 2, 0)

    rcx, rcy, rw, rh, rangle = mir
    rect = get_rect_corners(rcx, rcy, rw, rh, rangle)

    # -------------------------------------------------------------------------
    # Step 4: Determine unit depth
    # Cap depth so two unit rings + corridor fit within the rectangle's short side.
    # -------------------------------------------------------------------------
    depth = avg_depth
    min_dim = min(rw, rh)
    max_depth = (min_dim - corr_w - 8) / 2
    if depth > max_depth:
        depth = max(12, max_depth)

    door_w = 4         # Standard door/entry width (ft)
    cs = depth + door_w  # Corner set-back: how far a corner unit extends along an edge

    # -------------------------------------------------------------------------
    # Step 5: Set up edge direction vectors for the inscribed rectangle
    # dir_x / dir_y are the rotated axes; edge_dirs are the 4 outward edge directions.
    # -------------------------------------------------------------------------
    ca, sa = math.cos(rangle), math.sin(rangle)
    dir_x  = ( ca,  sa)
    dir_y  = (-sa,  ca)
    edge_dirs    = [dir_x, dir_y, vneg(dir_x), vneg(dir_y)]
    edge_normals = [dir_y, vneg(dir_x), vneg(dir_y), dir_x]  # Inward-pointing normals
    edge_lens    = [rw, rh, rw, rh]

    # -------------------------------------------------------------------------
    # Step 6: Place outer corner units (one at each rectangle corner)
    # Each corner unit spans `cs` along both edges meeting at that corner.
    # -------------------------------------------------------------------------
    corners = []
    for ci in range(4):
        v  = rect[ci]
        o0 = v
        o1 = vadd(v,  vscale(edge_dirs[ci], cs))
        i0 = vadd(o0, vscale(edge_normals[ci], depth))
        i1 = vadd(o1, vscale(edge_normals[ci], depth))
        corners.append(([o0, o1, i1, i0], 'corner'))

    # -------------------------------------------------------------------------
    # Step 7: Build the interleaved unit type queue
    # Units are distributed round-robin across types, scaled to units-per-floor.
    # -------------------------------------------------------------------------
    queue = []
    type_counts = [
        (u['type'], u['width_ft'], max(1, round(u['count'] / d_floors)))
        for u in data['dwelling_units']
    ]
    max_c = max(c for _, _, c in type_counts)
    for r in range(max_c):
        for utype, w, c in type_counts:
            if r < c:
                queue.append((utype, w))

    # -------------------------------------------------------------------------
    # Step 8: Place outer regular units along each edge
    # After subtracting the corner set-back, the remaining length is filled
    # with evenly-spaced units. Unit width is adjusted so they all fit cleanly.
    # -------------------------------------------------------------------------
    units = []
    qi = 0  # Queue index, wraps around
    units_by_edge = {}

    for ei in range(4):
        edge_start = len(units)
        avail = edge_lens[ei] - cs  # Available length after corner unit

        if avail < 10:
            units_by_edge[ei] = []
            continue

        origin = vadd(rect[ei], vscale(edge_dirs[ei], cs))
        uw = queue[qi % len(queue)][1] if queue else 20
        num = max(1, round(avail / uw))
        actual_w = avail / num

        # Widen units if they'd be too narrow (minimum 14 ft)
        while actual_w < 14 and num > 1:
            num -= 1
            actual_w = avail / num

        for k in range(num):
            o0 = vadd(origin, vscale(edge_dirs[ei],  k       * actual_w))
            o1 = vadd(origin, vscale(edge_dirs[ei], (k + 1)  * actual_w))
            i0 = vadd(o0, vscale(edge_normals[ei], depth))
            i1 = vadd(o1, vscale(edge_normals[ei], depth))
            utype = queue[qi % len(queue)][0]
            qi += 1
            units.append(([o0, o1, i1, i0], utype))

        units_by_edge[ei] = list(range(edge_start, len(units)))

    # -------------------------------------------------------------------------
    # Step 9: Compute core geometry
    # The "inner" rectangle (inside both unit rings) is the corridor rect.
    # Subtracting another corridor_width layer gives the service core.
    # -------------------------------------------------------------------------
    inner_w = rw - 2 * depth
    inner_h = rh - 2 * depth
    corr_rect = get_rect_corners(rcx, rcy, inner_w, inner_h, rangle) \
                if inner_w > 0 and inner_h > 0 else None

    core_w = inner_w - 2 * corr_w
    core_h = inner_h - 2 * corr_w
    raw_core = get_rect_corners(rcx, rcy, core_w, core_h, rangle) \
               if core_w > 0 and core_h > 0 else None

    # Chamfer the core corners for cleaner rendering
    if core_w > 0 and core_h > 0:
        core_chamfer = min(corr_w * 1.5, core_w / 4, core_h / 4)
        core_rect = []
        for ci2 in range(4):
            v   = raw_core[ci2]
            nxt = raw_core[(ci2 + 1) % 4]
            prv = raw_core[(ci2 - 1) % 4]
            dn = dist(v, nxt)
            dp = dist(v, prv)
            dn_dir = ((nxt[0] - v[0]) / dn, (nxt[1] - v[1]) / dn)
            dp_dir = ((prv[0] - v[0]) / dp, (prv[1] - v[1]) / dp)
            core_rect.append(vadd(v, vscale(dp_dir, core_chamfer)))
            core_rect.append(vadd(v, vscale(dn_dir, core_chamfer)))
    else:
        core_rect = None

    # -------------------------------------------------------------------------
    # Step 10: Determine inner layout mode based on core size
    # -------------------------------------------------------------------------
    inner_units    = []
    inner_corners  = []
    courtyard_rect = None
    boh_room       = None
    circ_room      = None
    layout_mode    = 'basic'
    inner_depth    = 0

    min_core = min(core_w, core_h) if (core_w > 0 and core_h > 0) else 0
    if min_core > 0:
        # Inner unit depth may be smaller than outer depth to leave room for a courtyard
        inner_depth = min(depth, (min_core - 2 * corr_w - 8) / 2)

    if inner_depth >= 12 and core_w > 0 and core_h > 0:
        # -----------------------------------------------------------------
        # Case 1: Courtyard — core is large enough for an inner unit ring
        # -----------------------------------------------------------------
        layout_mode = 'courtyard'
        inner_cs = min(inner_depth + door_w, min(core_w, core_h) / 2)

        # Inner corner units
        for ci in range(4):
            v  = raw_core[ci]
            o0 = v
            o1 = vadd(v,  vscale(edge_dirs[ci], inner_cs))
            i0 = vadd(o0, vscale(edge_normals[ci], inner_depth))
            i1 = vadd(o1, vscale(edge_normals[ci], inner_depth))
            inner_corners.append(([o0, o1, i1, i0], 'corner'))

        # Inner regular units (same logic as outer ring)
        inner_edge_lens    = [core_w, core_h, core_w, core_h]
        inner_units_by_edge = {}

        for ei in range(4):
            edge_start = len(inner_units)
            avail = inner_edge_lens[ei] - inner_cs

            if avail < 10:
                inner_units_by_edge[ei] = []
                continue

            origin = vadd(raw_core[ei], vscale(edge_dirs[ei], inner_cs))
            uw = queue[qi % len(queue)][1] if queue else 20
            num = max(1, round(avail / uw))
            actual_w = avail / num

            while actual_w < 14 and num > 1:
                num -= 1
                actual_w = avail / num

            for k in range(num):
                o0 = vadd(origin, vscale(edge_dirs[ei],  k       * actual_w))
                o1 = vadd(origin, vscale(edge_dirs[ei], (k + 1)  * actual_w))
                i0 = vadd(o0, vscale(edge_normals[ei], inner_depth))
                i1 = vadd(o1, vscale(edge_normals[ei], inner_depth))
                utype = queue[qi % len(queue)][0]
                qi += 1
                inner_units.append(([o0, o1, i1, i0], utype))

            inner_units_by_edge[ei] = list(range(edge_start, len(inner_units)))

        # Remaining open space after inner units = the courtyard
        cy_w = core_w - 2 * inner_depth - 2 * corr_w
        cy_h = core_h - 2 * inner_depth - 2 * corr_w
        courtyard_rect = get_rect_corners(rcx, rcy, cy_w, cy_h, rangle) \
                         if cy_w > 0 and cy_h > 0 else None

    elif core_w > 10 and core_h > 8:
        # -----------------------------------------------------------------
        # Case 2: core_rooms — BOH + circulation rooms placed in a full core
        # -----------------------------------------------------------------
        layout_mode = 'core_rooms'
        boh_w, boh_h   = 4, 4
        circ_w, circ_h = 6, 4
        total_w = boh_w + circ_w  # = 10

        if core_w >= core_h:
            # Place rooms side-by-side along the long axis
            boh_dx  = -total_w / 2 + boh_w  / 2
            circ_dx =  total_w / 2 - circ_w / 2
            boh_room  = get_rect_corners(rcx + boh_dx  * ca, rcy + boh_dx  * sa, boh_w,  boh_h,  rangle)
            circ_room = get_rect_corners(rcx + circ_dx * ca, rcy + circ_dx * sa, circ_w, circ_h, rangle)
        else:
            # Rotate rooms 90° to fit the taller core
            boh_dy  = -total_w / 2 + boh_w  / 2
            circ_dy =  total_w / 2 - circ_w / 2
            boh_room  = get_rect_corners(rcx - boh_dy  * sa, rcy + boh_dy  * ca, boh_h,  boh_w,  rangle)
            circ_room = get_rect_corners(rcx - circ_dy * sa, rcy + circ_dy * ca, circ_h, circ_w, rangle)

    elif core_w > 4 and core_h > 4:
        # -----------------------------------------------------------------
        # Case 3: notch — tight core; rooms are clipped to fit
        # -----------------------------------------------------------------
        layout_mode = 'notch'
        boh_w, boh_h   = 4, 4
        circ_w, circ_h = 6, 4
        total_w = boh_w + circ_w

        if core_w >= 10:
            # Fit side-by-side horizontally, clipping height to core
            boh_dx  = -total_w / 2 + boh_w  / 2
            circ_dx =  total_w / 2 - circ_w / 2
            boh_room  = get_rect_corners(rcx + boh_dx  * ca, rcy + boh_dx  * sa,
                                         boh_w,  min(boh_h,  core_h * 0.9), rangle)
            circ_room = get_rect_corners(rcx + circ_dx * ca, rcy + circ_dx * sa,
                                         circ_w, min(circ_h, core_h * 0.9), rangle)
        elif core_h >= 10:
            # Fit side-by-side vertically, clipping width to core
            boh_dy  = -total_w / 2 + boh_w  / 2
            circ_dy =  total_w / 2 - circ_w / 2
            boh_room  = get_rect_corners(rcx - boh_dy  * sa, rcy + boh_dy  * ca,
                                         min(boh_h,  core_w * 0.9), boh_w,  rangle)
            circ_room = get_rect_corners(rcx - circ_dy * sa, rcy + circ_dy * ca,
                                         min(circ_h, core_w * 0.9), circ_w, rangle)
        else:
            # Core is too small for two rooms — place a single small BOH room only
            boh_room = get_rect_corners(rcx, rcy,
                                        min(4, core_w * 0.9), min(4, core_h * 0.9), rangle)

    # -------------------------------------------------------------------------
    # Step 11: Fallback BOH + circulation room placement
    # If no rooms were placed above (basic mode, or courtyard without core rooms),
    # attach them to the inner edge of the most populated unit edges instead.
    # -------------------------------------------------------------------------
    if boh_room is None:
        boh_w, boh_h   = 4, 4
        circ_w, circ_h = 6, 4

        def _room_on_unit(unit_verts, rw, rh, on_outer):
            """
            Places a rectangular room on the outer or inner edge of a unit, centered.

            Args:
                unit_verts : (o0, o1, i1, i0) — the four corners of the host unit
                rw         : room width (along the edge)
                rh         : room depth (into the unit)
                on_outer   : if True, place on the street-facing edge; else corridor-facing

            Returns:
                List of 4 (x, y) corners for the room polygon.
            """
            o0, o1, i1, i0 = unit_verts
            if on_outer:
                p0, p1, into = o0, o1, i0
            else:
                p0, p1, into = i0, i1, o0

            e_len = dist(p0, p1)
            e_dir = ((p1[0] - p0[0]) / e_len, (p1[1] - p0[1]) / e_len)
            d_len = dist(p0, into)
            d_dir = ((into[0] - p0[0]) / d_len, (into[1] - p0[1]) / d_len)

            # Center the room along the edge
            off = max(0, (e_len - rw) / 2)
            s = vadd(p0, vscale(e_dir, off))

            return [
                s,
                vadd(s, vscale(e_dir, rw)),
                vadd(vadd(s, vscale(e_dir, rw)), vscale(d_dir, rh)),
                vadd(s, vscale(d_dir, rh))
            ]

        if layout_mode == 'courtyard' and inner_units:
            # Place on the outer (street-facing) edge of the two most-populated inner edges
            ranked = sorted(range(4), key=lambda e: len(inner_units_by_edge.get(e, [])), reverse=True)
            long_edges = [e for e in ranked if inner_units_by_edge.get(e)][:2]

            for idx, (rw2, rh2) in enumerate([(boh_w, boh_h), (circ_w, circ_h)]):
                ei  = long_edges[idx % len(long_edges)]
                ulist = inner_units_by_edge[ei]
                ui  = ulist[len(ulist) // 2]  # Middle unit on the edge
                room = _room_on_unit(inner_units[ui][0], rw2, rh2, on_outer=True)
                if idx == 0: boh_room  = room
                else:        circ_room = room

        elif units:
            # Place on the corridor-facing edge of the two most-populated outer edges
            ranked = sorted(range(4), key=lambda e: len(units_by_edge.get(e, [])), reverse=True)
            long_edges = [e for e in ranked if units_by_edge.get(e)][:2]

            for idx, (rw2, rh2) in enumerate([(boh_w, boh_h), (circ_w, circ_h)]):
                ei  = long_edges[idx % len(long_edges)]
                ulist = units_by_edge[ei]
                ui  = ulist[len(ulist) // 2]  # Middle unit on the edge
                room = _room_on_unit(units[ui][0], rw2, rh2, on_outer=False)
                if idx == 0: boh_room  = room
                else:        circ_room = room

    # -------------------------------------------------------------------------
    # Step 12: Return all computed geometry
    # -------------------------------------------------------------------------
    return {
        'boundary':       boundary,
        'rect':           rect,
        'corr_rect':      corr_rect,
        'core_rect':      core_rect,
        'units':          units,
        'corners':        corners,
        'depth':          depth,
        'corr_w':         corr_w,
        'fp_area':        fp_area,
        'pid':            pid,
        'rw':             rw,
        'rh':             rh,
        'avg_depth':      avg_depth,
        'layout_mode':    layout_mode,
        'inner_units':    inner_units,
        'inner_corners':  inner_corners,
        'courtyard_rect': courtyard_rect,
        'boh_room':       boh_room,
        'circ_room':      circ_room,
        'core_w':         core_w,
        'core_h':         core_h,
        'inner_depth':    inner_depth,
        'rcx':            rcx,
        'rcy':            rcy,
        'rangle':         rangle,
    }

'''
# â”€â”€ Text grid rasterizer â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
def rasterize(layout, cell_size=3):
    """Rasterize layout to text grid. Returns (grid_lines, diagnostics)."""
    bd = layout['boundary']
    allx=[p[0] for p in bd]; ally=[p[1] for p in bd]
    x0,x1 = min(allx)-5, max(allx)+5
    y0,y1 = min(ally)-5, max(ally)+5
    cols = int((x1-x0)/cell_size)+1
    rows = int((y1-y0)/cell_size)+1

    layout_mode = layout.get('layout_mode', 'basic')
    courtyard = layout.get('courtyard_rect')
    boh = layout.get('boh_room')
    circ = layout.get('circ_room')

    # Build spaces: corners first (priority), then regular, then inner
    spaces = []
    for verts,ctype in layout['corners']:
        spaces.append((verts, 'K', 'corner'))
    for verts,ctype in layout.get('inner_corners', []):
        spaces.append((verts, 'K', 'corner'))
    for verts,utype in layout['units']:
        spaces.append((verts, UNIT_CHAR.get(utype,'?'), utype))
    for verts,utype in layout.get('inner_units', []):
        spaces.append((verts, UNIT_CHAR.get(utype,'?'), utype))

    rect = layout['rect']
    corr = layout['corr_rect']
    core = layout['core_rect']

    grid = []
    unit_cells = {}
    corridor_cells = 0
    core_cells = 0
    courtyard_cells = 0
    gap_cells = 0
    overlap_cells = 0
    outside_cells = 0

    unit_corridor_adj = {}

    for r in range(rows):
        row = []
        y = y1 - r * cell_size
        for c in range(cols):
            x = x0 + c * cell_size

            in_rect = point_in_poly((x,y), rect)
            in_boundary = point_in_poly((x,y), bd)

            if not in_rect and not in_boundary:
                row.append(' ')
                outside_cells += 1
                continue

            if not in_rect and in_boundary:
                row.append('~')
                continue

            in_core = core and point_in_poly((x,y), core)
            in_corr = corr and point_in_poly((x,y), corr) and not in_core

            # BOH/Circ checked first â€” they overlay units
            if boh and point_in_poly((x,y), boh):
                row.append('B')
                continue
            if circ and point_in_poly((x,y), circ):
                row.append('C')
                continue

            claims = []
            for si,(verts,ch,label) in enumerate(spaces):
                if point_in_poly((x,y), verts):
                    claims.append((si, ch))

            if len(claims) > 1:
                corner_claims = [c for c in claims if c[1]=='K']
                if corner_claims:
                    si,ch = corner_claims[0]
                    row.append(ch)
                    unit_cells[ch] = unit_cells.get(ch,0)+1
                    if si not in unit_corridor_adj: unit_corridor_adj[si]=False
                    for dx2,dy2 in [(-1,0),(1,0),(0,-1),(0,1),(-1,-1),(1,-1),(-1,1),(1,1)]:
                        nx2,ny2=x+dx2*cell_size,y+dy2*cell_size
                        if corr and point_in_poly((nx2,ny2),corr) and not (core and point_in_poly((nx2,ny2),core)):
                            unit_corridor_adj[si]=True
                        if layout_mode=='courtyard' and core and point_in_poly((nx2,ny2),core):
                            if not (courtyard and point_in_poly((nx2,ny2),courtyard)):
                                unit_corridor_adj[si]=True
                else:
                    row.append('!')
                    overlap_cells += 1
            elif len(claims) == 1:
                si, ch = claims[0]
                row.append(ch)
                unit_cells[ch] = unit_cells.get(ch,0)+1
                if si not in unit_corridor_adj:
                    unit_corridor_adj[si] = False
                for dx2,dy2 in [(-1,0),(1,0),(0,-1),(0,1),(-1,-1),(1,-1),(-1,1),(1,1)]:
                    nx2,ny2 = x+dx2*cell_size, y+dy2*cell_size
                    if corr and point_in_poly((nx2,ny2),corr) and not (core and point_in_poly((nx2,ny2),core)):
                        unit_corridor_adj[si] = True
                    if layout_mode=='courtyard' and core and point_in_poly((nx2,ny2),core):
                        if not (courtyard and point_in_poly((nx2,ny2),courtyard)):
                            unit_corridor_adj[si] = True
            elif courtyard and point_in_poly((x,y), courtyard):
                row.append('O')
                courtyard_cells += 1
            elif boh and point_in_poly((x,y), boh):
                row.append('B')
            elif circ and point_in_poly((x,y), circ):
                row.append('C')
            elif in_core:
                if layout_mode == 'courtyard':
                    row.append('.')  # inner corridor
                    corridor_cells += 1
                else:
                    row.append('#')
                    core_cells += 1
            elif in_corr:
                row.append('.')
                corridor_cells += 1
            else:
                row.append('?')
                gap_cells += 1

        grid.append(''.join(row))

    # Diagnostics
    diag = []
    diag.append(f"Grid: {cols}x{rows} @ {cell_size}ft/cell | Mode: {layout_mode}")
    diag.append(f"Unit cells: {unit_cells}")
    diag.append(f"Corridor cells: {corridor_cells}, Core cells: {core_cells}")
    if courtyard_cells: diag.append(f"Courtyard cells: {courtyard_cells}")
    diag.append(f"Gaps: {gap_cells}, Overlaps: {overlap_cells}")

    no_access = []
    for si,(verts,ch,label) in enumerate(spaces):
        if si in unit_corridor_adj and not unit_corridor_adj[si]:
            no_access.append(f"  {label} unit #{si} ({ch}) â€” NO corridor access!")
    if no_access:
        diag.append("ACCESS ISSUES:")
        diag.extend(no_access)
    else:
        diag.append("All units have corridor access âœ“")

    return grid, diag
'''

def rasterize(layout, cell_size=3):
    """
    Converts a computed floor plan layout into a 2D text grid for visualization
    and diagnostics.

    Each cell in the grid represents a `cell_size` x `cell_size` ft square.
    The grid is populated by sampling a point at each cell's top-left corner
    and testing which space (if any) it falls inside.

    Character legend:
        ' '  - Outside the boundary entirely
        '~'  - Inside the boundary but outside the inscribed rectangle (irregular parcel area)
        'K'  - Corner unit
        '.'  - Corridor (outer ring or inner courtyard ring)
        '#'  - Service core (non-courtyard modes)
        'O'  - Open courtyard
        'B'  - Back-of-house (BOH) room
        'C'  - Circulation room
        '?'  - Gap (inside rect but unclaimed — indicates a layout issue)
        '!'  - Overlap (multiple non-corner units claimed the same cell)
        Letter codes for unit types (e.g. 'S', '1', '2', '3') via UNIT_CHAR lookup.

    Args:
        layout    : dict returned by compute_layout()
        cell_size : grid resolution in feet per cell (default: 3)

    Returns:
        grid  : list of strings, one per row, forming the text grid
        diag  : list of diagnostic message strings including cell counts
                and any units that lack corridor access
    """

    # -------------------------------------------------------------------------
    # Step 1: Compute grid dimensions from the boundary bounding box (+5 ft pad)
    # -------------------------------------------------------------------------
    bd = layout['boundary']
    all_x = [p[0] for p in bd]
    all_y = [p[1] for p in bd]
    x0, x1 = min(all_x) - 5, max(all_x) + 5
    y0, y1 = min(all_y) - 5, max(all_y) + 5
    cols = int((x1 - x0) / cell_size) + 1
    rows = int((y1 - y0) / cell_size) + 1

    # -------------------------------------------------------------------------
    # Step 2: Extract layout regions
    # -------------------------------------------------------------------------
    layout_mode = layout.get('layout_mode', 'basic')
    courtyard   = layout.get('courtyard_rect')
    boh         = layout.get('boh_room')
    circ        = layout.get('circ_room')

    rect = layout['rect']       # Inscribed rectangle corners
    corr = layout['corr_rect']  # Corridor rectangle (inside unit ring)
    core = layout['core_rect']  # Service core (inside corridor)

    # -------------------------------------------------------------------------
    # Step 3: Build the ordered list of occupiable spaces
    # Priority order matters: corners are tested before regular/inner units
    # so that overlapping corner cells resolve cleanly.
    # Each entry: (polygon_verts, display_char, label)
    # -------------------------------------------------------------------------
    spaces = []
    for verts, ctype in layout['corners']:
        spaces.append((verts, 'K', 'corner'))
    for verts, ctype in layout.get('inner_corners', []):
        spaces.append((verts, 'K', 'corner'))
    for verts, utype in layout['units']:
        spaces.append((verts, UNIT_CHAR.get(utype, '?'), utype))
    for verts, utype in layout.get('inner_units', []):
        spaces.append((verts, UNIT_CHAR.get(utype, '?'), utype))

    # -------------------------------------------------------------------------
    # Step 4: Initialize grid and diagnostic counters
    # -------------------------------------------------------------------------
    grid             = []
    unit_cells       = {}     # char -> cell count, for reporting unit coverage
    corridor_cells   = 0
    core_cells       = 0
    courtyard_cells  = 0
    gap_cells        = 0
    overlap_cells    = 0
    outside_cells    = 0

    # Tracks whether each space (by index into `spaces`) is adjacent to a corridor.
    # Used to flag units with no corridor access.
    unit_corridor_adj = {}

    # -------------------------------------------------------------------------
    # Step 5: Rasterize — iterate every cell and classify it
    # Grid origin is top-left; y decreases downward (row 0 = max y).
    # -------------------------------------------------------------------------
    for r in range(rows):
        row = []
        y = y1 - r * cell_size  # World-space Y (decreasing downward)

        for c in range(cols):
            x = x0 + c * cell_size  # World-space X

            in_rect     = point_in_poly((x, y), rect)
            in_boundary = point_in_poly((x, y), bd)

            # -- Cells fully outside all geometry --
            if not in_rect and not in_boundary:
                row.append(' ')
                outside_cells += 1
                continue

            # -- Irregular parcel area (between boundary and inscribed rect) --
            if not in_rect and in_boundary:
                row.append('~')
                continue

            # -- Determine which structural zones this cell is in --
            in_core = core and point_in_poly((x, y), core)
            in_corr = corr and point_in_poly((x, y), corr) and not in_core

            # -- BOH and Circ rooms overlay units — check them first --
            if boh  and point_in_poly((x, y), boh):
                row.append('B')
                continue
            if circ and point_in_poly((x, y), circ):
                row.append('C')
                continue

            # -- Find all spaces (units/corners) that claim this cell --
            claims = [
                (si, ch)
                for si, (verts, ch, label) in enumerate(spaces)
                if point_in_poly((x, y), verts)
            ]

            # -- Helper: check 8-directional neighbors for corridor adjacency --
            def check_corridor_adj(si):
                """Mark space si as corridor-adjacent if any neighbor cell is in the corridor."""
                if si not in unit_corridor_adj:
                    unit_corridor_adj[si] = False
                for dx, dy in [(-1,0),(1,0),(0,-1),(0,1),(-1,-1),(1,-1),(-1,1),(1,1)]:
                    nx, ny = x + dx * cell_size, y + dy * cell_size
                    # Outer corridor adjacency
                    if corr and point_in_poly((nx, ny), corr) and \
                       not (core and point_in_poly((nx, ny), core)):
                        unit_corridor_adj[si] = True
                    # Inner corridor adjacency (courtyard mode: core ring = inner corridor)
                    if layout_mode == 'courtyard' and core and point_in_poly((nx, ny), core):
                        if not (courtyard and point_in_poly((nx, ny), courtyard)):
                            unit_corridor_adj[si] = True

            if len(claims) > 1:
                # -- Overlap: prefer corner claims; otherwise flag as conflict --
                corner_claims = [c for c in claims if c[1] == 'K']
                if corner_claims:
                    si, ch = corner_claims[0]
                    row.append(ch)
                    unit_cells[ch] = unit_cells.get(ch, 0) + 1
                    check_corridor_adj(si)
                else:
                    row.append('!')
                    overlap_cells += 1

            elif len(claims) == 1:
                # -- Normal case: exactly one space claims this cell --
                si, ch = claims[0]
                row.append(ch)
                unit_cells[ch] = unit_cells.get(ch, 0) + 1
                check_corridor_adj(si)

            elif courtyard and point_in_poly((x, y), courtyard):
                row.append('O')
                courtyard_cells += 1

            elif in_core:
                if layout_mode == 'courtyard':
                    row.append('.')   # Inner courtyard corridor
                    corridor_cells += 1
                else:
                    row.append('#')   # Solid service core
                    core_cells += 1

            elif in_corr:
                row.append('.')
                corridor_cells += 1

            else:
                # Inside the rect but no space claimed it — layout gap
                row.append('?')
                gap_cells += 1

        grid.append(''.join(row))

    # -------------------------------------------------------------------------
    # Step 6: Build diagnostics
    # -------------------------------------------------------------------------
    diag = []
    diag.append(f"Grid: {cols}x{rows} @ {cell_size}ft/cell | Mode: {layout_mode}")
    diag.append(f"Unit cells: {unit_cells}")
    diag.append(f"Corridor cells: {corridor_cells}, Core cells: {core_cells}")
    if courtyard_cells:
        diag.append(f"Courtyard cells: {courtyard_cells}")
    diag.append(f"Gaps: {gap_cells}, Overlaps: {overlap_cells}")

    # Flag any units that were never found adjacent to a corridor
    no_access = [
        f"  {label} unit #{si} ({ch}) — NO corridor access!"
        for si, (verts, ch, label) in enumerate(spaces)
        if si in unit_corridor_adj and not unit_corridor_adj[si]
    ]
    if no_access:
        diag.append("ACCESS ISSUES:")
        diag.extend(no_access)
    else:
        diag.append("All units have corridor access ✓")

    return grid, diag

# â”€â”€ PNG renderer â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
def render_png(layout, ax):
    bd=layout['boundary'];rect=layout['rect']
    corr=layout['corr_rect'];core=layout['core_rect']
    courtyard=layout.get('courtyard_rect')
    boh=layout.get('boh_room')
    circ=layout.get('circ_room')
    layout_mode=layout.get('layout_mode','basic')

    ax.add_patch(MplPoly(bd,closed=True,fc='#fafafa',ec='#ccc',lw=0.8,zorder=-1))
    if corr: ax.add_patch(MplPoly(corr,closed=True,fc='#f0f0f0',ec='none',zorder=0))
    if core:
        if layout_mode == 'courtyard':
            ax.add_patch(MplPoly(core,closed=True,fc='#f0f0f0',ec='#999',lw=0.5,zorder=1))
        else:
            ax.add_patch(MplPoly(core,closed=True,fc='#d5d5d5',ec='#999',lw=0.5,zorder=1))
            cx,cy=centroid(core)
            ax.text(cx,cy,'CORE',ha='center',va='center',fontsize=7,fontweight='bold',color='#555',zorder=6)

    if courtyard:
        ax.add_patch(MplPoly(courtyard,closed=True,fc='#e8f5e9',ec='#4caf50',lw=1,zorder=2))
        cx,cy=centroid(courtyard)
        ax.text(cx,cy,'YARD',ha='center',va='center',fontsize=6,fontweight='bold',color='#2e7d32',zorder=6)

    if boh:
        ax.add_patch(MplPoly(boh,closed=True,fc='#ffccbc',ec='#333',lw=0.8,zorder=8))
        cx,cy=centroid(boh)
        ax.text(cx,cy,'BOH',ha='center',va='center',fontsize=5,fontweight='bold',color='#bf360c',zorder=9)
    if circ:
        ax.add_patch(MplPoly(circ,closed=True,fc='#d1c4e9',ec='#333',lw=0.8,zorder=8))
        cx,cy=centroid(circ)
        ax.text(cx,cy,'CIR',ha='center',va='center',fontsize=5,fontweight='bold',color='#4527a0',zorder=9)

    # Outer units
    for verts,ctype in layout['corners']:
        ax.add_patch(MplPoly(verts,closed=True,fc=COLORS['corner'],ec='#333',lw=0.5,zorder=2))
        ux=sum(v[0] for v in verts)/len(verts);uy=sum(v[1] for v in verts)/len(verts)
        ax.text(ux,uy,'CR',ha='center',va='center',fontsize=4,color='#888',zorder=3)
    for verts,utype in layout['units']:
        ax.add_patch(MplPoly(verts,closed=True,fc=COLORS.get(utype,'#eee'),ec='#333',lw=0.5,zorder=2))
        ux=sum(v[0] for v in verts)/len(verts);uy=sum(v[1] for v in verts)/len(verts)
        ax.text(ux,uy,LABELS.get(utype,'?'),ha='center',va='center',fontsize=5,color='#555',zorder=3)

    # Inner units (courtyard mode)
    for verts,ctype in layout.get('inner_corners', []):
        ax.add_patch(MplPoly(verts,closed=True,fc=COLORS['corner'],ec='#333',lw=0.5,zorder=4))
        ux=sum(v[0] for v in verts)/len(verts);uy=sum(v[1] for v in verts)/len(verts)
        ax.text(ux,uy,'CR',ha='center',va='center',fontsize=4,color='#888',zorder=5)
    for verts,utype in layout.get('inner_units', []):
        ax.add_patch(MplPoly(verts,closed=True,fc=COLORS.get(utype,'#eee'),ec='#333',lw=0.5,zorder=4))
        ux=sum(v[0] for v in verts)/len(verts);uy=sum(v[1] for v in verts)/len(verts)
        ax.text(ux,uy,LABELS.get(utype,'?'),ha='center',va='center',fontsize=5,color='#555',zorder=5)

    rx=[p[0] for p in rect]+[rect[0][0]];ry=[p[1] for p in rect]+[rect[0][1]]
    ax.plot(rx,ry,'k-',lw=2,zorder=5)
    if corr:
        cx2=[p[0] for p in corr]+[corr[0][0]];cy2=[p[1] for p in corr]+[corr[0][1]]
        ax.plot(cx2,cy2,color='#aaa',ls='--',lw=0.4,zorder=4)
    ax.set_aspect('equal');ax.axis('off')
    allx=[p[0] for p in bd];ally=[p[1] for p in bd];pad=20
    ax.set_xlim(min(allx)-pad,max(allx)+pad);ax.set_ylim(min(ally)-pad,max(ally)+pad)
    nu=len(layout['units'])+len(layout.get('inner_units',[]))
    nc=len(layout['corners'])+len(layout.get('inner_corners',[]))
    mode=layout.get('layout_mode','basic')
    ax.set_title(f"{layout['pid'].upper()} â€” {mode}\n{nu}+{nc}cr | {layout['rw']:.0f}x{layout['rh']:.0f} | d={layout['depth']:.0f}",fontsize=8,fontweight='bold')


# â”€â”€ SolverResult JSON output â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
TYPE_NAMES = {'studio':'Studio','1br':'1BR','2br':'2BR','3br':'3BR','corner':'Corner'}

def _make_space(sid, stype, name, floor_idx, verts, is_vert=False):
    area = poly_area(verts)
    return {
        'id':sid,'type':stype,'name':name,'floor_index':floor_idx,
        'geometry':{'vertices':[[round(v[0],2),round(v[1],2)] for v in verts]},
        'target_area_sf':round(area,1),'actual_area_sf':round(area,1),
        'membership':1,'area_deviation':'+0.0%','is_vertical':is_vert,
    }

def _make_rect_space(sid, stype, name, floor_idx, cx, cy, w, h, rot_deg=0, is_vert=False):
    area = w*h
    return {
        'id':sid,'type':stype,'name':name,'floor_index':floor_idx,
        'geometry':{'x':round(cx,2),'y':round(cy,2),'width':round(w,2),'height':round(h,2),'rotation':round(rot_deg,2)},
        'target_area_sf':round(area,1),'actual_area_sf':round(area,1),
        'membership':1,'area_deviation':'+0.0%','is_vertical':is_vert,
    }

def _local_to_global(rcx, rcy, ca, sa, lx, ly):
    return (rcx + lx*ca - ly*sa, rcy + lx*sa + ly*ca)

def _parking_spaces(data, layout, floor_idx):
    """Parking bays + support rooms."""
    spaces=[]
    rcx,rcy,rangle=layout['rcx'],layout['rcy'],layout['rangle']
    ca,sa=math.cos(rangle),math.sin(rangle)
    rot=math.degrees(rangle)
    hw,hh=layout['rw']/2,layout['rh']/2

    # Parking bays: 5 stalls wide (45ft) Ã— 2 deep (36ft)
    park=data.get('parking',{})
    total_stalls=park.get('underground_stalls',45)
    BAY_W,BAY_H,STALLS_PER=45,36,10
    n_bays=max(1,math.ceil(total_stalls/STALLS_PER))
    aisle_w=min(2*hh-30,80)

    placed=0;bc=0
    for row_sign in [1,-1]:
        row_y=row_sign*(aisle_w/2+BAY_H/2)
        x=-hw+10+BAY_W/2
        while x<hw-10 and placed<total_stalls:
            bc+=1
            stalls_here=min(STALLS_PER,total_stalls-placed)
            gx,gy=_local_to_global(rcx,rcy,ca,sa,x,row_y)
            spaces.append(_make_rect_space(f'parking_bay_{bc}_f{floor_idx}','PARKING',
                          f'P{placed+1}-{placed+stalls_here}',floor_idx,gx,gy,BAY_W,BAY_H,rot))
            placed+=stalls_here; x+=BAY_W+2
        if placed>=total_stalls: break

    # Support rooms â€” fit within the rect, wrapping rows if needed
    support=[('storage','SUPPORT','Storage',15,12),
             ('trash_recycle','SUPPORT','Trash/Recycle',12,10),
             ('fan_room','SUPPORT','Fan Room',12,10),
             ('fire_pump','SUPPORT','Fire Pump',10,10),
             ('mpoe','SUPPORT','MPOE',10,8)]
    margin=5
    sx=-hw+margin+support[0][3]/2
    sy=-hh+margin+support[0][4]/2
    for sid,stype,name,sw,sh in support:
        if sx+sw/2>hw-margin:  # wrap to next row
            sx=-hw+margin+sw/2
            sy+=max(s[4] for s in support)+2
        gx,gy=_local_to_global(rcx,rcy,ca,sa,sx,sy)
        spaces.append(_make_rect_space(f'{sid}_f{floor_idx}',stype,name,floor_idx,gx,gy,sw,sh,rot))
        sx+=sw+2
    return spaces

def _ground_spaces(data, layout, floor_idx):
    """Ground floor: lobby, leasing, amenities, support."""
    spaces=[]
    rcx,rcy,rangle=layout['rcx'],layout['rcy'],layout['rangle']
    ca,sa=math.cos(rangle),math.sin(rangle)
    rot=math.degrees(rangle)
    hw,hh=layout['rw']/2,layout['rh']/2
    s=min(1,min(hw,hh)/80)  # scale for small buildings

    rooms=[
        ('lobby','CIRCULATION','Lobby',       30*s,20*s,  0,       hh*0.4),
        ('leasing','SUPPORT','Leasing',       18*s,15*s,  hw*0.4,  hh*0.3),
        ('mail','SUPPORT','Mail/Package',     15*s,12*s, -hw*0.4,  hh*0.3),
        ('lounge','AMENITY','Lounge',         30*s,25*s, -hw*0.35, 0),
        ('fitness','AMENITY','Fitness',       25*s,20*s, -hw*0.25,-hh*0.35),
        ('restroom_m','SUPPORT','Restroom M', 12*s,10*s,  hw*0.4, -hh*0.2),
        ('restroom_f','SUPPORT','Restroom F', 12*s,10*s,  hw*0.4, -hh*0.35),
        ('trash','SUPPORT','Trash',           10*s,8*s,   hw*0.35, 0),
        ('bike','SUPPORT','Bike Storage',     20*s,18*s,  hw*0.25,-hh*0.45),
    ]
    for sid,stype,name,rw2,rh2,lx,ly in rooms:
        gx,gy=_local_to_global(rcx,rcy,ca,sa,lx,ly)
        spaces.append(_make_rect_space(f'{sid}_f{floor_idx}',stype,name,floor_idx,gx,gy,rw2,rh2,rot))
    return spaces

def _residential_spaces(data, layout, floor_idx):
    spaces=[];uid=0
    def add(verts,utype,pfx):
        nonlocal uid; uid+=1
        spaces.append(_make_space(f'{pfx}_{utype}_{uid}_f{floor_idx}','DWELLING_UNIT',
                                  f'{TYPE_NAMES.get(utype,utype)} {uid}',floor_idx,verts))
    # Skip corners â€” only regular units
    for v,t in layout['units']: add(v,t,'unit')
    for v,t in layout.get('inner_units',[]): add(v,t,'iunit')

    if layout['core_rect'] and layout.get('layout_mode')!='courtyard':
        spaces.append(_make_space(f'core_f{floor_idx}','CORE','Core',floor_idx,layout['core_rect']))
    if layout['boh_room']:
        spaces.append(_make_space(f'boh_f{floor_idx}','SUPPORT','BOH',floor_idx,layout['boh_room']))
    if layout['circ_room']:
        spaces.append(_make_space(f'circ_f{floor_idx}','CIRCULATION','Circ',floor_idx,layout['circ_room']))
    if layout.get('courtyard_rect'):
        spaces.append(_make_space(f'courtyard_f{floor_idx}','AMENITY','Courtyard',floor_idx,layout['courtyard_rect']))
    return spaces

def write_solver_json(data, layout, output_path):
    pid=data['project_id']
    fp_area=data['building']['floor_plate_sf']
    sd=data['building']['story_distribution']
    n_dwell=max(1,round(sd.get('dwelling',1)))
    n_mixed=max(1,round(sd.get('mixed_use',0)))
    n_park=max(0,round(sd.get('parking',0)))

    boundary=[[round(p[0],2),round(p[1],2)] for p in layout['boundary']]
    floors=[]

    for i in range(n_park):
        fi=-(n_park-i)
        floors.append({'floor_index':fi,'floor_type':'PARKING','boundary':boundary,
                       'area_sf':fp_area,'spaces':_parking_spaces(data,layout,fi)})
    for i in range(n_mixed):
        ft='GROUND' if i==0 else 'RESIDENTIAL_TYPICAL'
        floors.append({'floor_index':i,'floor_type':ft,'boundary':boundary,
                       'area_sf':fp_area,'spaces':_ground_spaces(data,layout,i)})
    for i in range(n_dwell):
        fi=n_mixed+i
        floors.append({'floor_index':fi,'floor_type':'RESIDENTIAL_TYPICAL','boundary':boundary,
                       'area_sf':fp_area,'spaces':_residential_spaces(data,layout,fi)})

    total=sum(len(f['spaces']) for f in floors)
    result={
        'success':True,'obstruction':0,'iterations':1,
        'message':f'Band layout v7 â€” {pid}','violations':[],
        'metrics':{'placement_rate':'100%','avg_membership':'1.0',
                   'total_spaces':total,'placed_spaces':total},
        'building':{'floors':floors,'stalks':[],
                    'metrics':{'total_floors':len(floors),'total_spaces':total,'cohomology_obstruction':0}},
    }
    with open(output_path,'w') as f: json.dump(result,f,indent=2)
    print(f"  JSON â†’ {output_path}")


# â”€â”€ Main â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
if __name__ == '__main__':
    import sys
    data_dir = Path(__file__).parent / 'web-viewer' / 'public' / 'data'
    out_dir = Path(__file__).parent / 'web-viewer' / 'output'
    out_dir.mkdir(exist_ok=True)
    do_json = '--json' in sys.argv

    fig,axes = plt.subplots(2,2,figsize=(14,14))

    for idx, pid in enumerate(['p1','p4','p7','p9']):
        with open(data_dir/f'{pid}_building.json') as f:
            data = json.load(f)

        print(f"\n{'='*60}")
        print(f"  {pid.upper()}: {data['project_name']}")
        print(f"{'='*60}")

        layout = compute_layout(data)
        print(f"  Rect: {layout['rw']:.0f}x{layout['rh']:.0f} ({layout['rw']*layout['rh']/layout['fp_area']*100:.0f}%)")
        print(f"  Depth: {layout['depth']:.0f}ft (avg from data: {layout['avg_depth']:.0f}ft)")
        print(f"  Core: {layout['core_w']:.0f}x{layout['core_h']:.0f} | Mode: {layout['layout_mode']}")
        print(f"  Outer: {len(layout['units'])} units + {len(layout['corners'])} corners")
        if layout['layout_mode'] == 'courtyard':
            print(f"  Inner: {len(layout['inner_units'])} units + {len(layout['inner_corners'])} corners (depth={layout['inner_depth']:.0f}ft)")

        if do_json:
            write_solver_json(data, layout, data_dir/f'{pid}_output.json')

        grid, diag = rasterize(layout, cell_size=3)
        if not do_json:
            print()
            for line in grid:
                print(f"  {line}")
        print()
        for d in diag:
            print(f"  {d}")

        render_png(layout, axes.flat[idx])

    fig.tight_layout(pad=2)
    out = out_dir/'band_layout_all.png'
    fig.savefig(out,dpi=100,bbox_inches='tight',facecolor='white')
    plt.close(fig)
    print(f"\nPNG: {out}")

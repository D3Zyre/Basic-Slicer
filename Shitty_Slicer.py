#from matplotlib import lines
import numpy
import math
from stl import mesh

def create_cube(dimension):
    "creates a dimension^3 cube with a corner at the origin and returns the mesh object"
    data = numpy.zeros(12, dtype=mesh.Mesh.dtype)

    # Top of the cube
    data['vectors'][0] = numpy.array([[0, 1, 1],
                                    [1, 0, 1],
                                    [0, 0, 1]])
    data['vectors'][1] = numpy.array([[1, 0, 1],
                                    [0, 1, 1],
                                    [1, 1, 1]])
    # Front face
    data['vectors'][2] = numpy.array([[1, 0, 0],
                                    [1, 0, 1],
                                    [1, 1, 0]])
    data['vectors'][3] = numpy.array([[1, 1, 1],
                                    [1, 0, 1],
                                    [1, 1, 0]])
    # Left face
    data['vectors'][4] = numpy.array([[0, 0, 0],
                                    [1, 0, 0],
                                    [1, 0, 1]])
    data['vectors'][5] = numpy.array([[0, 0, 0],
                                    [0, 0, 1],
                                    [1, 0, 1]])
    # Bottom of the cube
    data['vectors'][6] = numpy.array([[0, 1, 0],
                                    [1, 0, 0],
                                    [0, 0, 0]])
    data['vectors'][7] = numpy.array([[1, 0, 0],
                                    [0, 1, 0],
                                    [1, 1, 0]])
    # Back face
    data['vectors'][8] = numpy.array([[0, 0, 0],
                                    [0, 0, 1],
                                    [0, 1, 0]])
    data['vectors'][9] = numpy.array([[0, 1, 1],
                                    [0, 0, 1],
                                    [0, 1, 0]])
    # Right face
    data['vectors'][10] = numpy.array([[0, 1, 0],
                                    [1, 1, 0],
                                    [1, 1, 1]])
    data['vectors'][11] = numpy.array([[0, 1, 0],
                                    [0, 1, 1],
                                    [1, 1, 1]])

    data['vectors'] *= float(dimension)

    mesh_to_slice = mesh.Mesh(data.copy())

    return mesh_to_slice

def slice_mesh(input_mesh, layer_height = 0.2, line_width = 0.4, max_xyz = (235, 235, 255), number_of_walls = 3, print_speeds = (10, 25, 50, 50, 150), retraction_distance = 8):
    """
    print speeds = (first layer, outer wall, inner walls, infill, movement)
    """

    # to convert an STL to GCODE, we will follow these steps:
    # 1. We "slice" our computations so that we are in R2 instead of R3. Essentially, we are going to compute the same thing for a bunch of different z-heights, starting at "0" and going up by the layer height.
    # for each slice:
    #     3. We find which triangle intersect with the current plane we are working with (the slice).
    #     4. Then we find the lines where each of those triangles intersects our plane.
    #     5. We have a bunch of straight lines, but of course they can be grouped into lines that connect in a closed loop, so this is what we do now, we group lines together in closed loops of connecting lines.
    #     6. Since STLs use triangles to make shapes, we might have multiple triangles that intersected our plane in a straight line. We don't need to split that into multiple lines if they're all parallel, so this step is dedicated to reducing any unnecessary lines (same end result but potentially less seperate lines in our groups)
    #     7. Now we have lines grouped together, but we need to know where the inside of each group is. For us humans that's really easy just by looking at it. Not so simple for a computer. This function finds the inside of each closed line group, and returns that as a list of normal directions (which pairs up with the lines in the line groups)
    #     At this point we're ready to start making the toolpath.
    #     8. We generate a toolpath for the printer, which has its own steps:
    #          9. We generate the walls. The first wall is the outer wall, then we go inside of each line group by utilizing the normal direction that points inside of the loop, to make each inner wall
    #          TODO still working on this and just commenting to make more sense and figure out that issue later so yeah

    z_heights = [i/100 for i in range(0, int(max_xyz[2]*100+layer_height*100), int(layer_height*100))][1:]
    for z_height in z_heights:
        print("\rcurrent height: {:.2f}mm".format(z_height), end = "")
        triangles_on_plane = find_intersecting_triangles(input_mesh, z_height)
        if any(triangles_on_plane):
            intersection_lines = find_intersection_lines(input_mesh, z_height, triangles_on_plane)
            line_groups = group_lines(intersection_lines)
            line_normals = orient_inside_of_groups(line_groups)
            toolpath = generate_toolpath(line_groups, line_normals, line_width, number_of_walls, print_speeds, layer_height, max_xyz, retraction_distance)
            raise SystemExit()

def find_intersecting_triangles(input_mesh, z_height):
    """
    Finds the triangles of a mesh that intersect with the plane representing the layer we are currently slicing at
    returns a list of true false values corresponding to the triangles of the mesh
    """

    # if a triangle sits entirely on the plane it can be completely ignored as the triangles surrounding it will eventually work and that checks out
    # to check if a triangle intersects with the plane, since the plane is the XY plane, we just need to check if any of the points are on opposite sides of, or on the plane\

    triangles_on_plane = [False for i in range(len(input_mesh.v0))]

    for triangle in range(len(input_mesh.v0)):
        is_above_plane = ((input_mesh.v0[triangle][2] > z_height), (input_mesh.v1[triangle][2] > z_height), (input_mesh.v2[triangle][2] > z_height)) # if any of the points are above the plane they are 1
        is_below_plane = ((input_mesh.v0[triangle][2] < z_height), (input_mesh.v1[triangle][2] < z_height), (input_mesh.v2[triangle][2] < z_height)) # if any of the points are above the plane they are 1
        is_on_plane = ((input_mesh.v0[triangle][2] == z_height), (input_mesh.v1[triangle][2] == z_height), (input_mesh.v2[triangle][2] == z_height)) # if any of the points are above the plane they are 1

        if ( any(is_on_plane) or ( any(is_above_plane) and any(is_below_plane) ) ) and not all(is_on_plane): # checks if either a point is on the plane or if one point is above the plane and another is below, ignores if all are on the plane
            triangles_on_plane[triangle] = True

    return triangles_on_plane


def find_intersection_lines(input_mesh, z_height, triangles_on_plane):
    """
    Using the information found in find_intersecting_triangles, this function converts those intersecting triangles into lines along the plane, the lines of intersection of each triangle
    returns a list of coordinate pair tuples
    """

    # for each triangle in the mesh, if it intersects the plane, we need to find where it intersects.
    # to do this, we need to know which points are above the plane and which are below the plane, and then just draw a line and do some linear interpolation math stuff to find the point where those lines intersect the plane.
    # the line we're looking for is going to be the line between the two points we obtain

    point_1 = (0, 0, 0) # the location of the first point of the line
    point_1_set = False # whether or not the first point has been set yet
    point_2_set = False # whether or not the second point has been set yet
    point_2 = (0, 0, 0) # the location of the second point of the line

    list_of_coords = list()

    for triangle in range(len(input_mesh.v0)):
        point_1_set = False # reset to false
        point_2_set = False # reset to false
        point_1 = (0, 0, 0) # if not set then it should still be 0, this is to help with troubleshooting
        point_2 = (0, 0, 0) # same as above
        if triangles_on_plane[triangle]:
            point_0_above = bool(input_mesh.v0[triangle][2] > z_height) # each of these is just a True/False value of whether said point is above or below the plane, pretty self-explanatory with the variable names
            point_1_above = bool(input_mesh.v1[triangle][2] > z_height)
            point_2_above = bool(input_mesh.v2[triangle][2] > z_height)
            point_0_below = bool(input_mesh.v0[triangle][2] < z_height)
            point_1_below = bool(input_mesh.v1[triangle][2] < z_height)
            point_2_below = bool(input_mesh.v2[triangle][2] < z_height)

            # check if any of the points are on the plane:
            if input_mesh.v0[triangle][2] == z_height:
                point_1 = tuple(input_mesh.v0[triangle])
                point_1_set = True
            if input_mesh.v1[triangle][2] == z_height:
                if point_1_set:
                    point_2 = tuple(input_mesh.v1[triangle])
                    point_2_set = True
                else:
                    point_1_set = True
                    point_1 = tuple(input_mesh.v1[triangle])
            if input_mesh.v2[triangle][2] == z_height:
                if point_1_set:
                    point_2 = tuple(input_mesh.v2[triangle])
                    point_2_set = True
                else:
                    point_1_set = True
                    point_1 = tuple(input_mesh.v2[triangle])
            # check each possible case for combination of points above and below the plane, these are all copy/pasted and only slightly changed so I'm sure there's a better way to do it but this works right now so I'm not touching it
            if point_0_above and point_1_below:
                z_fraction = 1 - (input_mesh.v0[triangle][2] - input_mesh.v1[triangle][2] - z_height) / (input_mesh.v0[triangle][2] - input_mesh.v1[triangle][2]) # fraction is the number from 0 to 1 that represents where the plane is relative to the two points from the triangle, 0 being the lower of the points
                x = input_mesh.v0[triangle][0] + z_fraction*(input_mesh.v1[triangle][0] - input_mesh.v0[triangle][0])
                y = input_mesh.v0[triangle][1] + z_fraction*(input_mesh.v1[triangle][1] - input_mesh.v0[triangle][1])
                if point_1_set:
                    point_2 = tuple([x, y, z_height])
                    point_2_set = True
                else:
                    point_1_set = True
                    point_1 = tuple([x, y, z_height])
            if point_0_above and point_2_below:
                z_fraction = 1 - (input_mesh.v0[triangle][2] - input_mesh.v2[triangle][2] - z_height) / (input_mesh.v0[triangle][2] - input_mesh.v2[triangle][2]) # fraction is the number from 0 to 1 that represents where the plane is relative to the two points from the triangle, 0 being the lower of the points
                x = input_mesh.v0[triangle][0] + z_fraction*(input_mesh.v2[triangle][0] - input_mesh.v0[triangle][0])
                y = input_mesh.v0[triangle][1] + z_fraction*(input_mesh.v2[triangle][1] - input_mesh.v0[triangle][1])
                if point_1_set:
                    point_2 = tuple([x, y, z_height])
                    point_2_set = True
                else:
                    point_1_set = True
                    point_1 = tuple([x, y, z_height])
            if point_1_above and point_2_below:
                z_fraction = 1 - (input_mesh.v1[triangle][2] - input_mesh.v2[triangle][2] - z_height) / (input_mesh.v1[triangle][2] - input_mesh.v2[triangle][2]) # fraction is the number from 0 to 1 that represents where the plane is relative to the two points from the triangle, 0 being the lower of the points
                x = input_mesh.v1[triangle][0] + z_fraction*(input_mesh.v2[triangle][0] - input_mesh.v1[triangle][0])
                y = input_mesh.v1[triangle][1] + z_fraction*(input_mesh.v2[triangle][1] - input_mesh.v1[triangle][1])
                if point_1_set:
                    point_2 = tuple([x, y, z_height])
                    point_2_set = True
                else:
                    point_1_set = True
                    point_1 = tuple([x, y, z_height])
            if point_0_below and point_1_above:
                z_fraction = 1 - (input_mesh.v1[triangle][2] - input_mesh.v0[triangle][2] - z_height) / (input_mesh.v1[triangle][2] - input_mesh.v0[triangle][2]) # fraction is the number from 0 to 1 that represents where the plane is relative to the two points from the triangle, 0 being the lower of the points
                x = input_mesh.v1[triangle][0] + z_fraction*(input_mesh.v0[triangle][0] - input_mesh.v1[triangle][0])
                y = input_mesh.v1[triangle][1] + z_fraction*(input_mesh.v0[triangle][1] - input_mesh.v1[triangle][1])
                if point_1_set:
                    point_2 = tuple([x, y, z_height])
                    point_2_set = True
                else:
                    point_1_set = True
                    point_1 = tuple([x, y, z_height])
            if point_0_below and point_2_above:
                z_fraction = 1 - (input_mesh.v2[triangle][2] - input_mesh.v0[triangle][2] - z_height) / (input_mesh.v2[triangle][2] - input_mesh.v0[triangle][2]) # fraction is the number from 0 to 1 that represents where the plane is relative to the two points from the triangle, 0 being the lower of the points
                x = input_mesh.v2[triangle][0] + z_fraction*(input_mesh.v0[triangle][0] - input_mesh.v2[triangle][0])
                y = input_mesh.v2[triangle][1] + z_fraction*(input_mesh.v0[triangle][1] - input_mesh.v2[triangle][1])
                if point_1_set:
                    point_2 = tuple([x, y, z_height])
                    point_2_set = True
                else:
                    point_1_set = True
                    point_1 = tuple([x, y, z_height])
            if point_1_below and point_2_above:
                z_fraction = 1 - (input_mesh.v2[triangle][2] - input_mesh.v1[triangle][2] - z_height) / (input_mesh.v2[triangle][2] - input_mesh.v1[triangle][2]) # fraction is the number from 0 to 1 that represents where the plane is relative to the two points from the triangle, 0 being the lower of the points
                x = input_mesh.v2[triangle][0] + z_fraction*(input_mesh.v1[triangle][0] - input_mesh.v2[triangle][0])
                y = input_mesh.v2[triangle][1] + z_fraction*(input_mesh.v1[triangle][1] - input_mesh.v2[triangle][1])
                if point_1_set:
                    point_2 = tuple([x, y, z_height])
                    point_2_set = True
                else:
                    point_1_set = True
                    point_1 = tuple([x, y, z_height])
            
            if not point_1_set: # this should never happen if find_intersecting_triangles works properly
                print("ERROR IN FIND_INTERSECTION_LINES: point 1 wasn't set")
                print(input_mesh.v0[triangle], input_mesh.v1[triangle], input_mesh.v2[triangle])
                print(point_1, point_2)
            """
            if not point_2_set: # this will happen if a triangle only has a single point on the plane, and in that case we can ignore it, that's why this is commented out
                print("ERROR IN FIND_INTERSECTION_LINES: point 2 wasn't set")
                print(input_mesh.v0[triangle], input_mesh.v1[triangle], input_mesh.v2[triangle])
                print(point_1, point_2)
            """

            if point_1_set and point_2_set:
                list_of_coords.append(tuple([point_1, point_2]))

    # cleanup:

    # removing "line" entries that are the same point twice:
    for line in list_of_coords:
        if all([line[0][i] == line[1][i] for i in range(3)]):
            list_of_coords.remove(line)
    
    # removing any lines that appear multiple times:
    removing = True
    while removing: # while at least one line has been removed per iteration (since a line may appear even more than twice idk)
        removing = False
        for line in list_of_coords:
            if list_of_coords.count(line) > 1: # if there's more than one line that's the same
                list_of_coords.remove(line) # remove one of them
                removing = True


    return list_of_coords


def group_lines(intersection_lines):
    """
    takes in the lines from find_intersection_lines() and groups them together into sections of connected lines, as a list of lists of points. each point connects to the next, and the last point connects to the first
    """

    groups = list() # this is the list where each entry is a group, and those groups are an ordered list of points that make a line in segments
    groups.append([intersection_lines[0][0]]) # start with the first point of the first line
    next_point = intersection_lines[0][1] # set the next point we're looking for to be the other end of the line we just started
    cannot_be = 0 # the next point to be added can't be the one we just referenced otherwise we'll hit a dead end
    points_left = list() # list of all the points, and the points get removed as they are added to a group
    for line in range(len(intersection_lines)):
        for point in range(len(intersection_lines[line])):
            points_left.append(tuple([line, point]))
    current_group = 0 # we start on the 0th group
    found_during_loop = True # this will let us know if matches are still being found or if there are no more matches for our current group
    
    while len(points_left) > 0: # while not all points have been connected yet
        found_during_loop = False # reset this to: no we haven't found a match in the for loop
        for p in points_left:
            point = intersection_lines[p[0]][p[1]] # the current point we're looking at is not p, p is just an identifier. The point is actually point.
            if all([point[i] == next_point[i] for i in range(3)]) and p[0] != cannot_be: # if the point we're looking at matches the last point we added to the current group
                groups[current_group].append(point) # add this point to the group
                points_left.remove(p) # and remove it from points left
                swap_point = (-int(p[1])) + 1 # 1 becomes 0 and 0 becomes 1
                points_left.remove((p[0], swap_point)) # remove the next point we're looking for from the current line, because we don't need it any more (and we don't want our code to "find" the same line we were just on, and get stuck)
                next_point = intersection_lines[p[0]][swap_point] # set the next point we're looking for
                cannot_be = p[0] # to be honest not sure why this variable is needed but my code doesn't work without it so it's here
                found_during_loop = True # a match has been found in the for loop
        if not found_during_loop: # if no match was found that means we need to move to the next group
            current_group += 1 # move to the next group
            groups.append([intersection_lines[points_left[0][0]][points_left[0][1]]])
            cannot_be = points_left[0][0]

    for group in range(len(groups)): # just checking to make sure the first and last point of a group match
        if any([groups[group][0][i] != groups[group][-1][i] for i in range(3)]):
            print("ERROR IN GROUP_LINES: the first and last point of a group should match, but don't") # if they don't print an error
        else:
            groups[group].pop() # if they do we can remove the last entry because it isn't needed

    return remove_line_group_clutter(groups)


def remove_line_group_clutter(line_groups):
    """
    removes straight lines that are unnecessarily split into multiple straight lines
    """
    g = -1
    for group in line_groups:
        g += 1 # iterate group number
        remove = list() # this is the list of indexes we need to remove from the current group
        last_point = group[0] # set the initial last point to be the first one
        last_slope = (group[0][1] - group[-1][1] ) / ( group[0][0] - group[-1][0] ) # set the initial corresponding slope
        p = 0
        for point in group[1:]:
            p += 1 # iterate point number
            line = tuple([point[i] - last_point[i] for i in range(3)]) # get line vector from the current and last points
            slope = line[1]/line[0] # get the line's slope

            if last_slope == slope: # if the current and last slope are the same then the last point is bad
                remove.append(p-1) # so we want to remove the last point

            last_point = point # save the last point
            last_slope = slope # save the last slope
        
        remove.reverse() # we have to delete points last to first so that we don't mess up our indexes
        [line_groups[g].pop(i) for i in remove] # remove all the points that are unnecessary
    
    return line_groups


def orient_inside_of_groups(line_groups):
    """
    takes in the line groups from group_lines() and returns a respective list of normal unit vectors (same format as the output of group_lines())
    """
    # we need to know where the inside of the groups are to orient which way the walls thickness is going, and then which side the infill goes on
    # to do this we just need to have a normal vector for each line segment, pointing inside of the group
    # to do that, I'm going to pick one of the two possible 2D normals randomly for the first line segment of each group, then go along the line of the group, assigning that same direction to each line
    # once that's done, I'll check, for each line, if the normal vector hits another line. If it doesn't for any of the lines, then I reverse the direction of all the normals, and now it should work.

    groups_orientation = list()
    current_group = -1

    for group in line_groups:
        current_group += 1
        group.append(group[0]) # readding the first point to the end of the line
        groups_orientation.append(list())
        last_point = group[0] # set to the first point
        first_line = True
        for current_point in group[1:]:
            # first, for the first line in a group, we need to turn our initial direction vector into a normal unit vector
            # then, for each subsequent line, we need to calculate the new normal direction by rotating the last normal vector based on the direction of the line
            # basically we need to rotate the direction of the last normal vector until it matches the direction of the current normal vector, then we note the orientation from the rotated one and reverse it
            line_vector = tuple([current_point[i] - last_point[i] for i in  range(3)]) # calculating the direction vector of the line
            line_unit_vector = get_unit_vector(line_vector)
            normal_unit_vector = get_unit_normal_vector_R2(line_vector) # get unit normal vector
            if first_line: # not all of this is actually exclusive to first line, only references to initial_direction are
                groups_orientation[current_group].append(normal_unit_vector) # adding the unit normal vector to the list
                first_line = False # we're not on the first line anymore
            else:
                angle_between_lines = angle_between_vectors(line_unit_vector, last_line_unit_vector) # this gets the angle in the path between the current and last line
                cross = cross_product(last_line_unit_vector, line_unit_vector) # we need the direction of rotation as well, so we cross product the vectors, if the result is +z then the rotation was counterclockwise, otherwise clockwise
                normal_unit_vector = rotate_vector_R2(last_normal, angle_between_lines*cross[2]) # get the normal vector as a unit vector
                groups_orientation[current_group].append(normal_unit_vector) # adding the unit normal vector to the list

            last_line_unit_vector = line_unit_vector
            last_point = current_point # so that our code will know what the last point was
            last_normal = normal_unit_vector # same but with normal
        
    # now that we've made a normal vector for each line in the group that all match up, we need to check that they actually point inside as opposed to outside
    # if they point outside we just need to flip the vectors
    current_group = -1
    for group in line_groups:
        current_group += 1
        last_point = group[0] # set to the first point
        point_number = -1
        direction_of_first_line = tuple([group[1][i] - group[0][i] for i in range(3)]) # this is just the first line, but yeah same thing
        direction_of_first_line_unit = get_unit_vector(direction_of_first_line)
        total_signed_angle = 0
        for current_point in group[1:]:
            point_number += 1
            # for each group, we need to look at the direction the first line is going
            # and the direction of rotation of the group as we move along the line (it should go through 2 pi radians)
            # from those two data points we can determine which way the normal vector on the first line should be oriented
            # then apply than to all the lines in the group
            # to get the direction of rotation we'll have to add up all the signed angles for each line, I can't think of any other way around that
            # if the rotation is ccw then the normals should be ccw to the lines, and vice versa
            line_vector = tuple([current_point[i] - last_point[i] for i in range(3)]) # calculating the direction vector of the line
            line_unit_vector = get_unit_vector(line_vector)
            angle_between_lines = angle_between_vectors(line_unit_vector, last_line_unit_vector)
            cross = cross_product(last_line_unit_vector, line_unit_vector)
            total_signed_angle += angle_between_lines*cross[2] # this should be  + or - 2pi
            
            last_line_unit_vector = line_unit_vector
            last_point = current_point


        correct_normal_direction = all([rotate_vector_R2(direction_of_first_line_unit, math.pi/2 * total_signed_angle/abs(total_signed_angle))[i] >= groups_orientation[current_group][0][i] - 0.00000001 and rotate_vector_R2(direction_of_first_line_unit, math.pi/2 * total_signed_angle/abs(total_signed_angle))[i] <= groups_orientation[current_group][0][i] + 0.00000001 for i in range(3)])
        if not correct_normal_direction:
            # then we just need to flip the sign on all our normals for that group, simple as that
            for normal in range(len(groups_orientation[current_group])):
                groups_orientation[current_group][normal] = tuple([-i for i in groups_orientation[current_group][normal]])


    return groups_orientation


angle_between_vectors = lambda vector1, vector2 : math.acos(sum([(vector1[i]*vector2[i]) for i in range(len(vector1))]))
cross_product = lambda vector1, vector2 : tuple([vector1[1]*vector2[2] - vector1[2]*vector2[1], vector1[2]*vector2[0] - vector1[0]*vector2[2], vector1[0]*vector2[1] - vector1[1]*vector2[0]])
rotate_vector_R2 = lambda vector, radians : tuple([math.cos(radians)*vector[0] - math.sin(radians)*vector[1], math.sin(radians)*vector[0] + math.cos(radians)*vector[1], vector[2]])
get_length = lambda vector : math.sqrt(sum([i**2 for i in vector]))


def get_unit_vector(vector):
    length = get_length(vector)
    return tuple([i / length for i in vector])


def get_unit_normal_vector_R2(vector, swap = False):
    """
    if the input vector is in R3 that's ok, the output will be in R3 as well, but the last digit won't be used for this function, it's 0

    this function is really dumb and needs to be made better
    """
    try:
        if swap:
            initial_direction = (1, 0)
        else:
            initial_direction = (-1, 0)

        is_R3 = (len(vector) == 3)
        if is_R3:
            vector = vector[0:2]
        line_unit_vector = get_unit_vector(vector)
        top_dot_product = sum([(line_unit_vector[i]*initial_direction[i]) for i in range(2)]) # do the numerator dot product of the projection, the denominator is 1 cuz unit vectors
        proj = tuple([top_dot_product * (line_unit_vector[p]) for p in range(2)]) # projecting the initial direction vector onto the first line
        normal_vector = tuple([initial_direction[i] - proj[i] for i in range(2)]) # getting the actual normal vector
        if is_R3:
            normal_vector = (normal_vector[0], normal_vector[1], 0)
        length_of_normal = math.sqrt(sum([i**2 for i in normal_vector])) # getting the length of the normal vector we calculated
        if length_of_normal == 0:
            raise ZeroDivisionError # TODO, fix please, this is stupid
        normal_unit_vector = tuple([i / length_of_normal for i in normal_vector]) # calculating the unit vector of the normal vector we calculated
    except ZeroDivisionError:
        if swap:
            initial_direction = (0, 1)
        else:
            initial_direction = (0, -1)

        top_dot_product = sum([(line_unit_vector[i]*initial_direction[i]) for i in range(2)]) # do the numerator dot product of the projection, the denominator is 1 cuz unit vectors
        proj = tuple([top_dot_product * (line_unit_vector[p]) for p in range(2)]) # projecting the initial direction vector onto the first line
        normal_vector = tuple([initial_direction[i] - proj[i] for i in range(2)]) # getting the actual normal vector
        if is_R3:
            normal_vector = (normal_vector[0], normal_vector[1], 0)
        length_of_normal = math.sqrt(sum([i**2 for i in normal_vector])) # getting the length of the normal vector we calculated
        normal_unit_vector = tuple([i / length_of_normal for i in normal_vector]) # calculating the unit vector of the normal vector we calculated

    return normal_unit_vector


def generate_toolpath(line_groups, line_normals, line_width, number_of_walls, print_speeds, layer_height, max_xyz, retraction_distance): # call subfunctions
    wall_line_groups = generate_walls(line_groups, line_normals, line_width, number_of_walls)
    infill_line_groups = []

    gcode = convert_path_to_gcode(wall_line_groups, infill_line_groups, print_speeds, line_width, layer_height, max_xyz, retraction_distance)



def generate_walls(line_groups, line_normals, line_width, number_of_walls):
    """
    for each wall, we need to offset the line path into the object by line_width. the outer wall is offset by half of the line_width.
    we need to keep the new line path each time to make it easy to generate the next wall and so that it can be passed to the infill generator later
    
    wall_line_groups work as follows: wall_line_groups[group#][wall#][point#],
    the outer wall is the first wall, then it goes inside
    """
    # for each wall/iteration, we need to offset the line path
    # to do this, I'm going to very simply move the lines in their normal directions by the line_width,
    # and then I need to calculate where those lines now intersect, first for adjacent lines,
    # then also for all lines, in case the concentric wall generation ends up at another wall, in which case we would stop making new walls as that part of the object is full
    # in case walls have an angle greater than pi, we will actually need to extend the lines so that they intersect, (simple grade 10 math)
    # in the end, all adjacent lines need to intersect to create a new closed path for a wall

    # finding the intersection between two lines is simple, convert them to linear equations (grade 10 math) and then solve their intersection
    # then just extend or shorten them to that point, repeat for each pair of adjacent lines

    # TODO: what about bottom and top layers? or maybe that should be taken care of in infill, or in a seperate function?

    wall_line_groups = list() # a list of line groups, one list of groups for each wall
    wall_line_groups_orientation = list() # similarly, a list of the line's normals... TODO: check if this is even needed lol
    offset = 0

    # preallocate the lists
    wall_line_groups = [[list() for wall in range(round(number_of_walls))] for group in range(len(line_groups))]
    wall_line_groups_orientation = [[list() for wall in range(round(number_of_walls))] for group in range(len(line_groups))]

    outer_wall = True # the first wall we make is the outer wall and we treat it differently
    for wall in range(round(number_of_walls)): # I'm just going to ignore any decimals in walls
        current_group = -1
        if not outer_wall: # outer wall is treated slightly differently than inner wall
            offset += line_width # offset inwards by a line width
        else: # if we are on the outer wall
            offset += 0.5 * line_width # offset inwards by half a line width so the outer half is at the right place
            outer_wall = False # no longer on the outer wall exception
        for line_group in line_groups: # for each group
            current_group += 1
            current_line = -1
            #line_group.append(line_group[0]) # add first point to the end of the line ### for some reason this was bad? not sure why
            last_point = tuple([line_group[0][i] + line_normals[current_group][current_line][i] for i in range(3)]) # set to the first point
            # for the first line, we don't have a last line vector, which is needed for getting the last slope
            last_line_vector = tuple([last_point[i] - line_group[-2][i] for i in range(3)]) # to fix that we just need to calculate it manually
            last_line_unit_vector = get_unit_vector(last_line_vector) # convert to unit vector
            for current_point in line_group[1:]: # in the current group, for each point/line
                current_line += 1
                normal = line_normals[current_group][current_line] # update the current normal
                line_points = (tuple([current_point[i] + normal[i] * offset for i in range(3)]), last_point) # a tuple of the two points for the current line
                line_vector = tuple([current_point[i] - last_point[i] for i in range(3)]) # get the "slope" of our straight line
                line_unit_vector = get_unit_vector(line_vector) # get the unit vector of the slope
                # so we need to use the normal direction of the current and last line
                # and using that we offset those lines by line_width in their normal direction
                # then we calculate their intersection
                # line_points = tuple([tuple([point[i] + normal[i]*offset for i in range(3)]) for point in line_points]) # offset the current line by the wall thickness
                # the vector of the line stays the same because the offset line has the same direction and magnitude
                current_slope = line_unit_vector[1]/line_unit_vector[0] # the slope of the current line
                last_slope = last_line_unit_vector[1]/last_line_unit_vector[0] # the slope of the last line
                intersection_x = float( ( current_slope*line_points[0][0] - line_points[0][1] - last_slope*line_points[1][0] + line_points[1][1] ) / ( current_slope - last_slope ) ) # the x coordinate of the intersection between the current and last line after offsetting
                print(line_points)
                print("x")
                if not (abs(intersection_x) <= 1000000000000000): # if we have a problem calculating intersection_x due to infinite slope on one line or no change in slope between lines then we can assign it manually
                    print("x error")
                    if current_slope == last_slope: # if the slope is the same then the x point is the same as it was before (but offset)
                        intersection_x = line_points[1][0]
                    else: # this should be because one of the slopes was infinite, we need to know which one
                        if not (abs(current_slope) <= 1000000000000000): # if the current slope is infinite then the x point is any x point along that line
                            intersection_x = line_points[0][0]
                        else: # otherwise it's along the other line
                            intersection_x = line_points[1][0]
                intersection_y = float( ( current_slope*(intersection_x - line_points[0][0]) + line_points[0][1] ) ) # see intersection_x
                print("y")
                if not (abs(intersection_y) <= 1000000000000000): # see above
                    print("y error")
                    if current_slope == last_slope:
                        intersection_y = current_point[1]
                    else:
                        intersection_y = current_point[1]
                    intersection_y = line_points[0][1]
                print(intersection_x, intersection_y)
                intersection_point = tuple([intersection_x, intersection_y]) # the point where the current and last line intersect
                # new that we have the intersection of the two lines, we need to extend or shorten them to that point (moving the second point of the last line and the first of the current)
                line_points = tuple([intersection_point, line_points[1]]) # update the current point with the intersection
                #wall_line_groups[wall][current_group].append(line_points[1]) # add the previous (now completed) line to the list
                wall_line_groups[current_group][wall].append(line_points[1]) # add the previous (now completed) line to the list
                #wall_line_groups_orientation[wall][current_group].append(normal) # add the current normal to the list (this list is currently one entry ahead of the one above)
                wall_line_groups_orientation[current_group][wall].append(normal) # add the current normal to the list (this list is currently one entry ahead of the one above)
                    
                last_point = current_point # save the last point
                last_line_unit_vector = line_unit_vector # save the last unit vector

    print(wall_line_groups)
    return wall_line_groups

    # TODO: currently the function above doesn't check if a wall that it's creating is intersecting with a wall that has already been created (that isn't the adjacent one), fix that I guess


def generate_infill(): # needs to know where innermost wall is
    pass


def convert_path_to_gcode(wall_line_groups, infill_line_groups, print_speeds, line_width, layer_height, max_xyz, retraction_distance): # give warning if path outside of max xy
    """
    convention will be to retract at the end of something, not before something,
    convention will be to use absolute coordinates

    wall_line_groups work as follows: wall_line_groups[group#][wall#][point#]
    print speeds = (first layer, outer wall, inner walls, infill, movement)
    """

    gcode_commands_list = list() # this will quite simply be converted to text later, but essentially it's the same thing, just easier to work with this

    try: # make and use persistent variable
        E = convert_path_to_gcode.E # check if E exists, if it does, set E to it
    except AttributeError:
        convert_path_to_gcode.E = 0 # if it doesn't, create it
        E = 0

    for line_group in wall_line_groups: # we need to add a special case for the first point of each group!!!
        outer_wall = True
        first_point = True # the first point needs to have a movement without extrusion
        for wall in line_group:
            for point in wall:
                if first_point:
                    first_point = False # no longer on the first point
                    speed = print_speeds[4] # movement speed
                    X, Y, Z = point
                    gcode_commands_list.append("G0 F{:d} X{:.3f} Y{:.3f} Z{:.3f} ".format(speed*60, X, Y, Z))
                else:
                    speed = print_speeds[2-int(outer_wall)] # wall speeds
                    X, Y, Z = point
                    volume_per_mm = layer_height * line_width # cubic mm per mm (so yeah it's square mm)
                    volume_of_filament_per_mm = math.pi * (1.75/2)**2 # same as above
                    distance = get_length(tuple([point[i] - last_point[i] for i in range(3)])) # get distance travelled (mm)
                    E += volume_per_mm * distance / volume_of_filament_per_mm # volume / mm print * mm print / volume * mm filament == mm filament
                    gcode_commands_list.append("G1 F{:d} X{:.3f} Y{:.3f} Z{:.3f} E{:.5f}".format(speed*60, X, Y, Z, E))

                last_point = point
            outer_wall = False
    
    #[print(i) for i in gcode_commands_list]

def write_gcode_to_file(): # this needs to go after each layer's gcode has been generated and spliced together, we also add start and end gcode here
    pass


def main():
    cube_mesh = create_cube(10)
    slice_mesh(cube_mesh, 0.2, 0.4, (235, 235, 255), 3, (10, 25, 50, 50, 180), 8)

if __name__ == "__main__":
    main()


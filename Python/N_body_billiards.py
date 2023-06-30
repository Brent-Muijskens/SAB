from numpy import cos, sin
import matplotlib.pyplot as plt

def construct_boundary(points):
    # prepare the boundary
    N = len(points)
    lines = []

    for i in range(N):

        p1 = points[i-1]
        p2 = [points[i][0] - p1[0], points[i][1]-p1[1]]
        lines.append((p2, p1))

    return lines

def line_plotter(hits, style='b--', new=False):
    if new:
        plt.figure(figsize=(8,8))
        plt.xticks([],[]);
        plt.yticks([],[]);
    
    cooldown = 0
    for p1, p2 in zip(hits[:-1], hits[1:]):
        #if cooldown ==  0:
            #if abs(abs(p1[0]- p2[0]) -1) < 10**(-5):
            #    cooldown += 1
            #    continue
            #elif abs(abs(p1[1]- p2[1]) -1) < 10**(-5):
            #    cooldown += 1
            #    continue
        #else:
        #    cooldown = 0
        plt.plot([p1[0], p2[0]], [p1[1], p2[1]], style, lw=3)

def collision(pos, line, wall):
        # determine the collision of two vectors
        det = (wall[0][1]*line[0] - wall[0][0]*line[1])
        if det == 0:
            return None, None

        vec = [pos[0]-wall[1][0], pos[1]-wall[1][1]]

        result = [1/det*(-line[1]*vec[0] + line[0]*vec[1]), 1/det*(-wall[0][1]*vec[0] +  wall[0][0]*vec[1])]

        if result[0] <= 1 and result[0] > 0 and result[1] > 0:
            
            #return [result[0]*wall[0][0]+wall[1][0], result[0]*wall[0][1]+wall[1][1]]
            return [result[1]*line[0]+pos[0], result[1]*line[1]+pos[1]], result[0]

        return None, None

def reflect(line, wall):
    """
    Calculates the reflection of a line in a wall
    """
    normal = [-wall[0][1], wall[0][0]]
    factor = 2/(normal[1]**2 + normal[0]**2)
    dot_prod = factor*(line[0]*normal[0]+line[1]*normal[1])
    r = [line[k] - dot_prod*normal[k] for k in range(2)]

    return r

def simulate(domain_points, N_particles, positions, angles, time_limit = 10, interactions = []):

    # simulation data containers
    walls = construct_boundary(domain_points)
    unconf_walls = [-1]*N_particles
    unconf_wall_times = [[] for _ in range(N_particles)]
    domino = [[] for _ in range(N_particles)]

    particle_times = [0]*N_particles

    # create vectors
    vectors = [[cos(angle), sin(angle)] for angle in angles]

    # reflect wall ind
    reflect_indices =  [None]*N_particles

    # keep track of all positions
    pos_hist = [[pos] for pos in positions]

    COUNTER = len(domain_points)
    while COUNTER < 1000:
        # choose particle with smallest current time
        current_time = min(particle_times)
        # set the active particle
        ap = particle_times.index(current_time)

        if current_time >= time_limit:
            #print(len(walls))
            break
    
        # update vector
        if reflect_indices[ap] != None:
            r = reflect(vectors[ap], walls[reflect_indices[ap]])
            if abs(r[0]+vectors[ap][0]) < 10**(-8) and abs(r[1]+vectors[ap][1]) < 10**(-8):
            #   print('perp hit')
                particle_times[ap] += time_limit
                continue
            else:
                vectors[ap] = r

        # set position/vector
        pos, vec = positions[ap], vectors[ap]

        # simulation
        hit_norm = 1000
        hit_position = []

        # reset unconf wall data
        unconf_walls[ap] = -1
        domino[ap] = []

        for ind, wall in enumerate(walls):

            # calculate new hit candidate
            hit_candidate, scaler = collision(pos, vec, wall)

            # check if valid hit
            if hit_candidate == None:
                continue

            # check if the collision is closer than a potential previous found collision
            cnorm = ((pos[0]-hit_candidate[0])**2 + (pos[1]-hit_candidate[1])**2)**0.5

            # hit is too close -> to avoid collisions with self
            if cnorm < 0.00001:
                continue
            
            # see if hit is closer than previous found hits
            if cnorm < hit_norm:
                hit_norm = cnorm
                hit_position, hit_ind, wall_scaler  = hit_candidate, ind, scaler

        # stop when hit is too close -> avoid floating point precision incidents
        if hit_norm < 0.0001:
            particle_times[ap] += time_limit
            continue
        
        # no hit found -> shouldn't happen more for debugging 
        if hit_position == []:
            particle_times[ap] += time_limit
            continue

        # add wall index for reflection
        reflect_indices[ap] = hit_ind
        

        # test if collision is with an `unconfirmed` wall
        # if yes see which particle was earlier and adjust the lines for it
        try:
            hit_part = unconf_walls.index(hit_ind)
            t_0, delta_t = unconf_wall_times[hit_part]

            collision_time = t_0 + wall_scaler*delta_t
            # if line should not be there
            if collision_time > particle_times[ap] + hit_norm:

                # update wall
                prev_vec = walls[hit_ind][0]
                walls[hit_ind][0] = [wall_scaler*prev_vec[0], wall_scaler*prev_vec[1]] 

                # reset reflect indices
                reflect_indices[ap] = None
                reflect_indices[hit_part] = COUNTER 

                # revert time/pos
                particle_times[hit_part] = collision_time
                positions[hit_part] = hit_position
                pos_hist[hit_part][-1] = hit_position
                unconf_wall_times[hit_part] = [t_0, wall_scaler*delta_t]
                domino[ap].append((hit_part, 1))

                # domino
                for ki, si in domino[hit_part]:
                    if si < wall_scaler:
                        continue
                    reflect_indices[ki] = None        
                domino[hit_part] = []
            else:
                domino[hit_part].append((ap, wall_scaler))

        except ValueError:
            pass

        # update position/vector
        positions[ap] = hit_position

        # store position for plots
        pos_hist[ap].append(hit_position)

        # add wall
        walls.append([[hit_position[0]-pos[0], hit_position[1]-pos[1]], pos])
        unconf_walls[ap] = COUNTER
        unconf_wall_times[ap] = [particle_times[ap], hit_norm]

        # update particle time    
        particle_times[ap] += hit_norm

        # increment counter
        COUNTER += 1

    return pos_hist
        

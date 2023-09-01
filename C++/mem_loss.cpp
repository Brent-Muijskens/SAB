#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <array>
#include <tuple>
#include <algorithm>
#include <ctime>
#include <sstream>
#include <chrono>

struct Wall{
    std::array<float, 2> vect;
    std::array<float, 2> offset;
    int id = -1;
    float t0 = 0;
    float te = 0;
};

struct Domino{
    std::vector<float> scalars;
    Domino(int N_particles){
        scalars.assign(N_particles, -1);
    }
};

struct Hit_data{
    std::array<float, 2> vect;
    float scalar = -1;
};

struct Particle_data{
    // indices to check
    std::vector<int> seq;
    // prev indices
    std::array<int,2> prev_ind;

    // wall orientations
    std::vector<int> orientations;
    // wall neigbours
    std::vector<std::array<int, 2>> neigbours;
};

std::vector<Wall> construct_boundary(std::vector<std::array<float, 2>> &points, int N){
    /* DECLARATIONS*/
    std::vector<Wall> result;

    std::array<float, 2> p1; std::array<float, 2> p2;
    Wall temp;

    for (int i = 0; i < N; i++){

        if (i == 0){
            p1 = points[N-1];
            p2[0] = points[i][0] - p1[0];
            p2[1] = points[i][1] - p1[1];
        }
        else{
            p1 = points[i-1];
            p2[0] = points[i][0] - p1[0];
            p2[1] = points[i][1] - p1[1];
        }
        temp.vect = p2; temp.offset = p1;
        result.push_back(temp);
    }
    return result;
}

std::vector<std::array<float,2>> sample_polygon(std::vector<std::array<float, 2>> &points, 
                                                std::vector<Wall> &lines, int N, int N_samples){
    /* Sample polygon with accept reject method*/

    /* Declarations*/
    std::vector<std::array<float, 2>> samples;

    std::array<float, 2> vec; std::array<float, 2> factors;

    float x_sample; float y_sample;
    /* hard coded but such large figures should not be used*/
    float x_min = 10000000; float x_max = -10000000;;
    float y_min = 10000000;; float y_max = -10000000;;

    float det;

    int N_hits = 0;
    int intersections = 0;

    /* find bouding rectangle*/
    for (int i = 0; i < N; i++){
        float temp_x = points[i][0];
        float temp_y = points[i][1];

        if (temp_x < x_min){
            x_min = temp_x;
        }
        else if (temp_x > x_max){
            x_max = temp_x;
        }

        if (temp_y < y_min){
            y_min = temp_y;
        }
        else if (temp_y > y_max){
            y_max = temp_y;
        }
    }

    while (N_hits < N_samples){
        while (intersections % 2 == 0){
            /* New sample */
            x_sample = ((float)rand() / (float)RAND_MAX)*(x_max-x_min) + x_min;
            y_sample = ((float)rand() / (float)RAND_MAX)*(y_max-y_min) + y_min;
            
            for (int i = 0; i < N; i++){
                det = -lines[i].vect[0];
                if (det == 0){
                    continue;
                }
                vec = {x_sample-lines[i].offset[0], y_sample-lines[i].offset[1]};
                factors = {-1/det*(vec[0]), 1/det*(-lines[i].vect[1]*vec[0] + lines[i].vect[0]*vec[1])};
                if (factors[0] <= 1 && factors[0] > 0 && factors[1] > 0 ){
                    intersections ++;
                }
            }
        }

        /* Add sample*/
        samples.push_back({x_sample, y_sample});
        /* increment number of samples found */
        N_hits ++;
        /* reset number of intersections */
        intersections = 0;
    }
    return samples;
}

std::array<float, 2>& reflect(std::array<float,2> &vec, Wall &wall){
    /* Calculate reflection of vector of a wall*/
    static std::array<float, 2> result;

    float temp = 2/(std::pow(wall.vect[0], 2) + std::pow(-wall.vect[1], 2));
    temp = temp*(-vec[0]*wall.vect[1] + vec[1]*wall.vect[0]);

    result = {vec[0] + temp*wall.vect[1], vec[1] - temp*wall.vect[0]};

    return result;
}

Hit_data& collision(std::array<float,2> &pos, std::array<float,2> &vec, Wall &wall, float det){
    float offset0 = pos[0] - wall.offset[0];
    float offset1 = pos[1] - wall.offset[1];

    std::array<float, 2> temp;
    static Hit_data result;
    
    temp[0] = 1/det*(-vec[1]*offset0 + vec[0]*offset1);
    temp[1] = 1/det*(-wall.vect[1]*offset0 + wall.vect[0]*offset1);

    if (temp[0] <= 1 && temp[0] > 0 && temp[1] > 0){
        result.vect = {temp[1]*vec[0] + pos[0], temp[1]*vec[1] + pos[1]};
        result.scalar = temp[0];

        return result;
    }
    result.scalar = -1;
    return result;
}

void simulate(std::vector<std::array<float, 2>> &points, 
                int N, int N_particles,
                std::vector<Wall> &init_walls, 
                std::vector<std::array<float, 2>> start_positions,
                std::vector<std::array<float, 2>> start_vectors,
                std::vector<std::vector<int>> interactions,
                std::ofstream &DataFile,
                std::ofstream &DistData,
                float decay_delay=10,
                int time_limit = 150){

    /* DECLARATIONS*/
  
    /* CONSTANTS*/
    static const float NORM_TOL = pow(10,-7); 
    static const float NORM_TOL2 = pow(10,-4);
    static const float BASE_NORM = pow(10,5);
    static const float PERP_ERR = pow(10,-8);
    static const int MAX_ITERS = 100000; 
    static const int TIME_LIMIT = time_limit;

    
    /* CONTAINERS */
    std::vector<int> unconf_walls(N_particles, -1);
    std::vector<std::array<float, 2>> unconf_walls_times(N_particles, {0,0});
 
    std::vector<float> particle_times(N_particles, 0);
    std::vector<float> Particle_Dist(N_particles, 0);
    std::vector<int> reflect_indices(N_particles, -1);

    static Domino reset_domino(N_particles);
    static std::vector<Domino> domino (N_particles, reset_domino);

    /* positions of particles */
    std::vector<std::array<float, 2>> positions = start_positions;
    std::vector<std::array<float, 2>> vectors = start_vectors;

    /* Data */
    //std::vector<std::vector<std::array<float, 2>>> POS_HIST;

    /* init the storage list*/
    //for (int i = 0; i < N_particles; i++){
    //    POS_HIST.push_back({start_positions[i]});
    //}

    /* init the boundary*/
    std::vector<Wall> walls = init_walls; 
    Wall wall;
    Wall new_wall;

    /* Misc variables*/
    float cnorm; 
    float hit_norm;
    std::array<float, 2> hit_position;
    int hit_ind; float wall_scalar;
    float t_0; float delta_t; float collision_time;
    std::array<float, 2> prev_vec;

    float current_time;
    int ap;     /*active particle*/

    /* position and vector variables*/
    std::array<float,2> pos;
    std::array<float,2> vec;

    
    int COUNTER = N;

    /*MAIN SIMULATION LOOP*/
    while (COUNTER < MAX_ITERS){
        /*Choose particle with smallest current time*/
        auto it = std::min_element(std::begin(particle_times), std::end(particle_times));
        int ap = std::distance(std::begin(particle_times), it);
        float current_time = *it;

        /* Time limit check */
        if (current_time >= TIME_LIMIT){
            break;
        }

        /* update vector */
        if (reflect_indices[ap] != -1){
            std::array<float, 2> r = reflect(vectors[ap], walls[reflect_indices[ap]]);
            /* Test for perpendicular hit with some error margin*/
            if (std::abs(r[0]+vectors[ap][0]) < PERP_ERR && std::abs(r[1]+vectors[ap][1]) < PERP_ERR){
                particle_times[ap] += TIME_LIMIT;
                continue;
            }
            else{
                vectors[ap] = r;
            }   
        }
        /* Set position/vector for current particle */
        pos = positions[ap]; vec = vectors[ap];

        /* Prepare for next simulation loop*/
        hit_norm = BASE_NORM;
        unconf_walls[ap] = -1;
     
        domino[ap] = reset_domino;

        for (int i = 0; i < COUNTER; i++){
            wall = walls[i];

            if (interactions[ap][wall.id] == 1 && unconf_walls[wall.id] != i) {
                continue;
            }

            /* Test if wall is still `active' */
            if (wall.id != -1 && wall.te + decay_delay < current_time && unconf_walls[wall.id] != i){
                continue;
            }

            float det = (wall.vect[1]*vec[0] - wall.vect[0]*vec[1]);
            if (det == 0){
                continue;
            }
            /* get hit location and corresponding scalar*/
            Hit_data hit_data = collision(pos, vec, wall, det);

            /* non valid hit found*/
            if (hit_data.scalar == -1){
                continue;
            }

            /* calculate norm of next hit*/
            cnorm = sqrt(pow(pos[0]-hit_data.vect[0], 2) + pow(pos[1]-hit_data.vect[1], 2));

            if (wall.id != -1 && unconf_walls[wall.id] != i && wall.t0 + hit_data.scalar*(wall.te-wall.t0) + decay_delay < current_time + cnorm){
                continue;
            }

            /* hit too close -> continue*/
            if (cnorm < NORM_TOL){
                continue;
            }

            if (cnorm < hit_norm){
                hit_norm = cnorm;
                hit_position = hit_data.vect; hit_ind = i; wall_scalar = hit_data.scalar;
            }
        }
        /* NO collisions found*/
        if (hit_norm == BASE_NORM){
            particle_times[ap] += TIME_LIMIT;
            continue;
        }
        /* add index for reflection */
        reflect_indices[ap] = hit_ind;
        
        /* test if collision is with unconfirmed wall*/
        for (int j = 0; j < N_particles; j++){
            if (hit_ind == unconf_walls[j]){
                t_0 = unconf_walls_times[j][0]; delta_t = unconf_walls_times[j][1];

                collision_time = t_0 + wall_scalar*delta_t;
                /* check if wall should not be there */
                if (collision_time > current_time + hit_norm && interactions[j][ap] == 0){
                    /*update wall*/
                    prev_vec = walls[hit_ind].vect;
                    walls[hit_ind].vect = {wall_scalar*prev_vec[0], wall_scalar*prev_vec[1]};
                    walls[hit_ind].te = collision_time;

                    /* reset reflect indices*/
                    reflect_indices[ap] = -1;
                    reflect_indices[j] = COUNTER;

                    /* revert time/pos */
                    particle_times[j] = collision_time;
                    positions[j] = hit_position;
                    unconf_walls_times[j] = {t_0, wall_scalar*delta_t};
                    //POS_HIST[j].back() = hit_position;
                    Particle_Dist[j] = collision_time;
                    domino[ap].scalars[j] = 1;

                    /* domino */
                    for (int k = 0; k < N_particles; k++){
                        if (domino[j].scalars[k] == -1){
                            continue;
                        }
                        else if (domino[j].scalars[k] < wall_scalar){
                            continue;
                        }
                        reflect_indices[k] = -1;
                    }
                    domino[j] = reset_domino;
                }
                else if (interactions[ap][j] == 1){
                    reflect_indices[ap] = -1;
                }
                else if (collision_time + decay_delay < current_time + hit_norm){
                    reflect_indices[ap] = -1;
                }
                else{
                    domino[j].scalars[ap] = wall_scalar;
                }
            }
        }

        /* update position*/
        positions[ap] = hit_position;
       
        /* add the new wall */
        new_wall.vect = {hit_position[0]-pos[0], hit_position[1]-pos[1]};
        new_wall.offset = pos;
        new_wall.id = ap;
        new_wall.t0 = current_time;
        new_wall.te = current_time + hit_norm;
        walls.push_back(new_wall);
        
        unconf_walls[ap] = COUNTER;
        unconf_walls_times[ap] = {current_time, hit_norm};
        /*update particle time*/
        particle_times[ap] += hit_norm;
        Particle_Dist[ap] += hit_norm;

        /* Store position */
        //POS_HIST[ap].push_back(hit_position);

        /* hit is too close -> continue to avoid floating point precision incidents*/
        if (hit_norm < NORM_TOL2 && reflect_indices[ap] != -1){
            particle_times[ap] += TIME_LIMIT;
        }

        /* increment counter */
        COUNTER ++;
    }

    /* SAVE DATA*/
    for (int i = 0; i < N_particles; i++){
        //int data_size = POS_HIST[i].size() -1 ;
            // Position data
            DataFile << positions[i][0] << "," << positions[i][1] << "; ";
            // Distance data
            DistData << Particle_Dist[i]  << ";";
        }
    DataFile << "\n";
    DistData << "\n";
}

int main()
{    
    // Start measuring time
    auto begin = std::chrono::high_resolution_clock::now();

    /* init random seed */
    srand(time(NULL)); rand();
    
    /* Define pi*/
    const float pi = 3.14159274101257324219;

    /* Data files */
    static std::ofstream DataFile;
    static std::ofstream DistData;
    DataFile.open("posdata.txt", std::ios::out);
    DistData.open("distdata.txt", std::ios::out);

    int N_particles;
    int N_SIMULATIONS;
    float decay_delay;
    std::vector<std::array<float, 2>> points;
    std::vector<std::vector<int>> interactions;

    /* Get settings from file */
    std::ifstream settings ("model_settings.txt");
    std::string line;
    if (settings.is_open()){ 
        while(std::getline(settings, line) ){
            if (line[0] == '#' || line.empty()){
                continue;
            }
            // position of = sign
            auto delimiterPos = line.find("=");
            auto name = line.substr(0, delimiterPos);
            auto value = line.substr(delimiterPos + 1);

            /* Assign values */
            if (name == "N_particles"){
                // store sim data
                DataFile << "# " << line << "\n";
                DistData << "# " << line << "\n";

                N_particles = std::stoi(value);
            
            }
            else if (name == "N_simulations"){
                // store sim data
                DataFile << "# " << line << "\n";
                DistData << "# " << line << "\n";

                N_SIMULATIONS = std::stoi(value);
                
            }
            else if (name == "Domain"){
                // store sim data
                DataFile << "# " << line << "\n";
                DistData << "# " << line << "\n";
                
                value.erase(remove(value.begin(), value.end(), '('), value.end());
                value.erase(remove(value.begin(), value.end(), ')'), value.end());
                value.erase(remove(value.begin(), value.end(), ','), value.end());

                std::istringstream iss(value);
                std::string val;
                int c = 0;
                float temp;
                while(std::getline(iss, val, ' ')) {
                    if (c == 0){
                        temp = std::stof(val);
                        c += 1;
                    }
                    else if (c == 1){
                        points.push_back({temp, std::stof(val)});
                        c = 0;
                    }
                }
            }
            else if (name == "Interactions"){

                // init the interaction vector
                std::vector<int> base_interactions(N_particles, 0);
                std::vector<std::vector<int>> temp_interactions(N_particles, base_interactions);

                // store sim data
                DataFile << "# " << line << "\n";
                DistData << "# " << line << "\n";

                value.erase(remove(value.begin(), value.end(), '('), value.end());
                value.erase(remove(value.begin(), value.end(), ','), value.end());
                value.erase(remove(value.begin(), value.end(), ')'), value.end());
                
                std::istringstream iss(value);
                std::string val;
                int count = 0;
                while(std::getline(iss, val, ' ')) {
                    if (val.empty()){
                        continue;
                    }
                    if (val != ";"){
                        temp_interactions[count][std::stoi(val)] = 1;
                        continue;
                    }
                    count ++;
                }
                interactions = temp_interactions;
            }
            else if (name == "Decay_delay"){
                DataFile << "# " << line << "\n";
                DistData << "# " << line << "\n";

                decay_delay = std::stof(value);
            }
        }
    }

    else{
        std::cout << "Could not find settings file" << std::endl;
        return 0;
    }
    
    // Setup vectors for position and direction data
    std::vector<std::array<float, 2>> start_positions;
    std::vector<std::array<float, 2>> start_vectors(N_particles, {0,0});
    
    // Stuff for status prints
    int N = points.size();
    int sim_status_data = N_SIMULATIONS/10;
    if (sim_status_data == 0){
        sim_status_data = 1;
    }

    std::vector<Wall> BOUNDARY_LINES = construct_boundary(points, N);
    
    std::cout << "Starting: " << N_SIMULATIONS << " simulations with " << N_particles << " particles " << std::endl << "Progress: ";
    for (int sim = 0; sim < N_SIMULATIONS; sim++){
        /* random start positions*/
        start_positions = sample_polygon(points, BOUNDARY_LINES, N, N_particles);
        /* random start vectors */
        for (int np = 0; np < N_particles; np++){
            float rand_angle = ((float)rand()/(float)RAND_MAX)*2*pi;
            start_vectors[np] = {std::cos(rand_angle), std::sin(rand_angle)};
        }

        if (sim % sim_status_data == 0){
            std::cout << 10*(sim/sim_status_data) << "%, ";
        }
        simulate(points, N, N_particles, BOUNDARY_LINES, start_positions, start_vectors, interactions, DataFile, DistData, decay_delay);
    }
    // status print
    std::cout << "100%" << std::endl;
    // close files
    DataFile.close();
    DistData.close();

    // Stop measuring time and calculate the elapsed time
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);

    std::cout << "Simulation completed after: " << elapsed.count() * 1e-9 << " seconds" << std::endl;
    return 0;
}
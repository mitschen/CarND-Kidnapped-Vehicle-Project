/*
 * particle_filter.h
 *
 * 2D particle filter class.
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#ifndef PARTICLE_FILTER_H_
#define PARTICLE_FILTER_H_

#include <random>
#include "helper_functions.h"

/**
 * observation entry of each particle
 *
 * An observation information is transformed in an SObservation object
 * for each particle. This observation object contains after different
 * particle operations
 * * the id of the landmark with best match
 * * the x-position of the observation in map coord-sys
 * * the y-position of the observation in map coord-sys
 * * the resulting distance to the landmark best match
 * * the square of distance in x between landmark and observation
 * * the square of distance in y between landmark and observation
 * The last two values are used for error/ weight calculation
 */
struct SObservation
{
  /**
   * default initializer necessary for using it in vector
   * @Please note: should never been called
   */
  SObservation()
  : id(-1), x(0.), y(0.), distanceTo(-1.), x_diff_quad(-1.), y_diff_quad(-1.){};
  /**
   * initialization used during observation step
   */
  SObservation(double const &mapX, double const &mapY)
  : id(-1), x(mapX), y(mapY), distanceTo(-1.), x_diff_quad(-1.), y_diff_quad(-1.){};
  int id;
  double x;
  double y;
  double distanceTo;  //resulting distance to the landmark
  double x_diff_quad; //diff between landmark in x, quadric
  double y_diff_quad; //diff between landmark in y, quadric
};


/**
 * a certain particle object
 *
 * A particle is initialized with a certain state information containing
 * x, y and theta (=yaw). There are different functions which allows to
 * proceed the particle object through the different staps of an particle
 * filter. Please refer to the function description below
 */
struct Particle {
  /**
   * default C'tor necessary to allow usage in std::vector
   */
  Particle();
  /**
   * initialize a particle with a certain id and a state vector
   * The state vector represents x, y and yaw (theta)
   */
  Particle(int const &_id, double const state[]);
  /**
   * predict new particle position based on dynamics and noise (particle filter step)
   *
   * This method will update the internal state of an particle according to
   * the veloctiy and yawrate, taking a noise information in x, y, theta
   * passed as std_pos into consideration.
   *
   * The function following the math provided in lesson 14 part 7.
   *
   * Please note: it is necessary to pass an external random-engine to each
   * of the particles. Otherwise, the gaussian noise calculation will result
   * in very strange values,
   */
  void predictPosition(double const &dt, double const std_pos[], double const &velocity, double const &yawrate, std::default_random_engine &gen);
  /**
   * transform the passed observation into map-representation according to
   * orientation and location of a certain particle (particle filter step)
   *
   * The function following the math provided in lesson 14 part 14.
   */
  void transformAndUpdateObservation(std::vector<LandmarkObs> const &);
  /**
   * calculate the closest neighbor of the particles observations in reference
   * to the real map data. Consider the maxDistance to filter out
   * map-information we don't need to take into considerations.
   */
  void matchObservationsToMap(Map const &map_landmarks, double const &maxDistance);
  /**
   * calculate the particles weight by multiplying all distance "errors"
   * of all observations. Store the weight in the corresponding member.
   * This method will furthermore take an external noise into consideration
   *
   * The function following the math provided in lesson 14 part 17.
   */
  void calculateWeight(double const std_landmark[]);
  /**
   * print the details of the particle to std::out
   */
  void printMe();
  //MEMBERS
	int id;       //id of the particle
	double x;     //current x position
	double y;     //current y position
	double theta; //current yaw-angle
	double weight;//resulting weight (is valid after calculateWeight function is called)
	//list of SObservations objects. List will be upated after each transformAndUpdateObservation
	std::vector<SObservation> observations;
};



class ParticleFilter {
	
	// Number of particles to draw
	int num_particles; 
	
	
	
	// Flag, if filter is initialized
	bool is_initialized;
	
	// Vector of weights of all particles
	std::vector<double> weights;
	
	//Generator we're using
	//It took me while to figure out that using local random-engine
	//for each noise calculation results in always the same results
	//That leads to a very strange behaviour of the particles.
	std::default_random_engine randomGen;
public:
	
	// Set of current particles
	std::vector<Particle> particles;

	// Constructor
	// @param num_particles Number of particles
	ParticleFilter(int noParticals = 5)
	: num_particles(noParticals)
	, is_initialized(false)
	, weights(0)
	, randomGen()
	, particles(0){}

	// Destructor
	~ParticleFilter() {}

	/**
	 * init Initializes particle filter by initializing particles to Gaussian
	 *   distribution around first position and all the weights to 1.
	 * @param x Initial x position [m] (simulated estimate from GPS)
	 * @param y Initial y position [m]
	 * @param theta Initial orientation [rad]
	 * @param std[] Array of dimension 3 [standard deviation of x [m], standard deviation of y [m]
	 *   standard deviation of yaw [rad]]
	 */
	void init(double x, double y, double theta, double std[]);

	/**
	 * prediction Predicts the state for the next time step
	 *   using the process model.
	 * @param delta_t Time between time step t and t+1 in measurements [s]
	 * @param std_pos[] Array of dimension 3 [standard deviation of x [m], standard deviation of y [m]
	 *   standard deviation of yaw [rad]]
	 * @param velocity Velocity of car from t to t+1 [m/s]
	 * @param yaw_rate Yaw rate of car from t to t+1 [rad/s]
	 */
	void prediction(double delta_t, double std_pos[], double velocity, double yaw_rate);
	
	/**
	 * dataAssociation Finds which observations correspond to which landmarks (likely by using
	 *   a nearest-neighbors data association).
	 * @param predicted Vector of predicted landmark observations
	 * @param observations Vector of landmark observations
	 */
	void dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations);
	
	/**
	 * updateWeights Updates the weights for each particle based on the likelihood of the 
	 *   observed measurements. 
	 * @param sensor_range Range [m] of sensor
	 * @param std_landmark[] Array of dimension 2 [Landmark measurement uncertainty [x [m], y [m]]]
	 * @param observations Vector of landmark observations
	 * @param map Map class containing map landmarks
	 */
	void updateWeights(double sensor_range, double std_landmark[], const std::vector<LandmarkObs> &observations,
			const Map &map_landmarks);
	
	/**
	 * resample Resamples from the updated set of particles to form
	 *   the new set of particles.
	 */
	void resample();

	/*
	 * Set a particles list of associations, along with the associations calculated world x,y coordinates
	 * This can be a very useful debugging tool to make sure transformations are correct and assocations correctly connected
	 */
	Particle SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y);
	
	std::string getAssociations(Particle best);
	std::string getSenseX(Particle best);
	std::string getSenseY(Particle best);

	/**
	 * initialized Returns whether particle filter is initialized yet or not.
	 */
	const bool initialized() const {
		return is_initialized;
	}
};



#endif /* PARTICLE_FILTER_H_ */

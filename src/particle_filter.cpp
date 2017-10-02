/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <sstream>
#include <string>
#include <iterator>
#include <cassert>

#include "particle_filter.h"

using namespace std;

Particle::Particle()
: id(-1), x(0.), y(0.), theta(0.), weight(1.), observations(){};//, associations(), sense_x(), sense_y(){}
//Particle::Particle(int const &_id, double const &_x, double const& _y, double const& _theta)
//: id(_id), x(_x), y(_y), theta(_theta), weight(1.), associations(), sense_x(), sense_y(){}
Particle::Particle(int const &_id,double const state[])
: id(_id), x(state[0]), y(state[1]), theta(state[2]), weight(1.), observations(){};//, associations(), sense_x(), sense_y(){}

void Particle::predictPosition(double const &dt, double const std_pos[], double const &velocity, double const &yawrate, default_random_engine &gen)
{
  if(yawrate == 0.)
  {
    double const vMt(velocity*dt);
    x = x + vMt*cos(theta);
    y = y + vMt*sin(theta);
    //theta = theta; - no change in theta
  }
  else
  {
    double const vDy(velocity / yawrate); //velocity Divided yawrate
    double const yMt(yawrate*dt); //yawrate Multiplied time
    //PLEASE NOTE: update theta first - this has influence
    //on x and y. If not doing so, i've always running into a
    //high error rate.
    theta = theta + yMt;
    x = x + vDy * (sin(theta + yMt) - sin(theta));
    y = y + vDy * (-cos(theta + yMt) + cos(theta));
  }

  //now adding noise
  double state[] = {x, y, theta};
  enrichWithRandomNoise(state, std_pos, gen);
  x = state[0];
  y = state[1];
  theta = state[2];
}

void Particle::transformAndUpdateObservation(std::vector<LandmarkObs> const &obs)
{
  //identify some constants for the current iteration
  double const cost(cos(theta));
  double const sint(sin(theta));
  observations.clear();
  //transform each observation into a map representation and store this
  //observation with the particle
  for(int i(0), _maxI(obs.size()); i < _maxI; i++)
  {
    double const &xObservation(obs[i].x);
    double const &yObservation(obs[i].y);
    double mapX(x + cost * xObservation - sint * yObservation);
    double mapY(y + sint * xObservation + cost * yObservation);
    observations.push_back(SObservation(mapX, mapY));
  }
}

void Particle::matchObservationsToMap(Map const &map_landmarks, double const &maxDistance)
{
  std::vector<Map::single_landmark_s> resultingCandidates;
  double const lowerX(x-maxDistance), higherX(x+maxDistance),
               lowerY(y-maxDistance), higherY(y+maxDistance);

  //filter all landmarks which could not be in range according to
  //particle position and sensor range
  for(int i(0), _maxI(map_landmarks.landmark_list.size()); i < _maxI; i++)
  {
    Map::single_landmark_s const &lm(map_landmarks.landmark_list[i]);
    if(   (lm.x_f >= lowerX)
        &&(lm.x_f <= higherX)
        &&(lm.y_f >= lowerY)
        &&(lm.y_f <= higherY)
        )
    {
      resultingCandidates.push_back(lm);
    }
  }

  //iterate over each landmark and each particle. Identify for
  //each particle the best matching landmark
  for(int i(0), _maxI(resultingCandidates.size()); i < _maxI; i++)
  {
    Map::single_landmark_s const &lm(resultingCandidates[i]);
    for(int j(0), _maxJ(observations.size()); j < _maxJ; j++)
    {
      SObservation &obs(observations[j]);
      double x_diff(obs.x-lm.x_f);
      double y_diff(obs.y-lm.y_f);
      double x_diff_quad(x_diff*x_diff);
      double y_diff_quad(y_diff*y_diff);
      double const distance(sqrt(x_diff_quad+y_diff_quad));
      if(obs.id==-1)
      {
        obs.id = lm.id_i;
        obs.distanceTo = distance;
        obs.x_diff_quad = x_diff_quad;
        obs.y_diff_quad = y_diff_quad;
      }
      else
      {
        if(obs.distanceTo>distance)
        {
          obs.id = lm.id_i;
          obs.distanceTo = distance;
          obs.x_diff_quad = x_diff_quad;
          obs.y_diff_quad = y_diff_quad;
        }
      }
    }
  }

  //in case that there could not be a landmark in range
  //set the weight to 0.
  if(resultingCandidates.size() == 0)
  {
    observations.clear();
    weight = 0.;
  }
}

void Particle::calculateWeight(double const std_landmark[])
{
  double const nominator(2*M_PI*std_landmark[0]*std_landmark[1]);
  double const nom2std_x(2*std_landmark[0]*std_landmark[0]);
  double const nom2std_y(2*std_landmark[1]*std_landmark[1]);

  //if there wasn't a single observation for the current particle,
  //we set the weight to zero
  if(observations.size()!=0)
  {
    weight = 1.;
  }
  else
  {
    weight = 0.;
  }
  for(int i(0), _maxI(observations.size()); i<_maxI; i++)
  {
    SObservation const &obs(observations[i]);
    weight = weight * (exp(- (obs.x_diff_quad/nom2std_x + obs.y_diff_quad/nom2std_y)) / nominator);
  }
}

void Particle::printMe()
{
  printf("Particle %i, x=%lf, y=%lf, theta=%lf, weight=%lf, noObs=%i\n", id, x, y, theta, weight, (int)observations.size());
  for(int i(0); i<observations.size(); i++)
  {
    SObservation const &obs(observations[i]);
    printf("Observation %i xdiff %lf ydiff %lf\n", obs.id, sqrt(obs.x_diff_quad), sqrt(obs.y_diff_quad));
  }
}

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).


  //following the example from Lesson14 / 4&5

  //some constants for the standard deviation using for initialization
//  double const &std_x(std[0]);
//  double const &std_y(std[1]);
//  double const &std_theta(std[2]);
//
//  //initialize the random generator
//  normal_distribution<double> dist_x(x, std_x);
//  normal_distribution<double> dist_y(y, std_y);
//  normal_distribution<double> dist_theta(theta, std_theta);
////  default_random_engine gen;

  double state[3];
  state[0]=x;
  state[1]=y;
  state[2]=theta;

  for(int i(0); i<num_particles; i++)
  {
//    state[0] = dist_x(randomGen);
//    state[1] = dist_y(randomGen);
//    state[2] = dist_theta(randomGen);
    enrichWithRandomNoise(state, std, randomGen);
    particles.push_back(Particle(i, state));
//    particles.push_back(Particle(i, dist_x(gen), dist_y(gen), dist_theta(gen)));
  }

  is_initialized=true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
  assert(initialized());
  for(int i(0); i < num_particles; i++)
  {
    particles[i].predictPosition(delta_t, std_pos, velocity, yaw_rate, randomGen);
  }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

  //Not used at all. I've followed the oop approach and realized the different steps
  //in the subclasses. See the updateWeights function
  assert(false);

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html


  //MScharf
  //as discussed in Section 14
  //Must be normalized between 0 and 1 because it is used as propabilities
  double overall_weight(0.);
  weights.clear();
  static int counter(1);
  for(int i(0); i<num_particles; i++)
  {
    Particle &particle(particles[i]);
    //associate the observation with each particle
    particle.transformAndUpdateObservation(observations);
    //match the particle observation the the "real" map landmarks
    particle.matchObservationsToMap(map_landmarks, sensor_range);
    //calculate the resulting weights (based on the distance-errors)
    particle.calculateWeight(std_landmark);
    //store the weights of the particles
    weights.push_back(particle.weight);
    overall_weight+=particle.weight;
  }
  std::cout << __FUNCTION__ <<" Iteration "<<counter<<" and overall weight of "<<overall_weight<<std::endl;
  counter++;
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
  std::random_device rd;
  std::mt19937 gen(rd());
  std::discrete_distribution<> d(weights.begin(), weights.end());
  std::vector<Particle> weighted_particles(0);
  for(int n=0; n<num_particles; ++n) {
    weighted_particles.push_back(particles[d(gen)]);
  }
  particles = weighted_particles;
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
  //I've followed the oop appraoch and pushed the different operations
  //to subclases. See the header for further information
  assert(false);
 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
  vector<int> v;
  //Workaround: don't touch the main-function.
  for(int i(0), _maxI(best.observations.size()); i < _maxI; i++)
  {
    v.push_back(best.observations[i].id);
  }
  stringstream ss;
  copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}
string ParticleFilter::getSenseX(Particle best)
{
  vector<int> v;
  //Workaround: don't touch the main-function.
  for(int i(0), _maxI(best.observations.size()); i < _maxI; i++)
  {
    v.push_back(best.observations[i].x);
  }
  stringstream ss;
  copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}
string ParticleFilter::getSenseY(Particle best)
{
  vector<int> v;
  //Workaround: don't touch the main-function.
  for(int i(0), _maxI(best.observations.size()); i < _maxI; i++)
  {
    v.push_back(best.observations[i].y);
  }
  stringstream ss;
  copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

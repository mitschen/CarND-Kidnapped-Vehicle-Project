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

void Particle::predictPosition(double const &dt, double const std_pos[], double const &velocity, double const &yawrate)
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
    //TODO: do we need to normalize the theta value/ yawrate multiplied with time?
    //      to a value between 0-2PI?
    x = x + vDy * (sin(theta + yMt) - sin(theta));
    y = y + vDy * (-cos(theta + yMt) + cos(theta));
    theta = theta + yMt;
  }

  //now adding noise
  double state[] = {x, y, theta};
  enrichWithRandomNoise(state, std_pos);
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
  for(int i(0), _maxI(obs.size()); i < _maxI; i++)
  {
    double const &xObservation(obs[i].x);
    double const &yObservation(obs[i].y);
    double mapX(x + cost * xObservation - sint * yObservation);
    double mapY(y + sint* + xObservation + cost * yObservation);
    observations.push_back(SObservation(mapX, mapY));
  }
}

void Particle::matchObservationsToMap(Map const &map_landmarks, double const &maxDistance)
{
  std::cout<<__FUNCTION__<<" starting with initial set of landmarks of size "<<map_landmarks.landmark_list.size()<<std::endl;
  std::vector<Map::single_landmark_s> resultingCandidates;
  double const lowerX(x-maxDistance), higherX(x+maxDistance),
               lowerY(y-maxDistance), higherY(y+maxDistance);
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
  std::cout<<__FUNCTION__<<" reducing landmarks to size "<<resultingCandidates.size()<<" according to sensor range "<<maxDistance<<std::endl;
  for(int i(0), _maxI(resultingCandidates.size()); i < _maxI; i++)
  {
    Map::single_landmark_s const &lm(resultingCandidates[i]);
    for(int j(0), _maxJ(observations.size()); j < _maxJ; j++)
    {
      SObservation &obs(observations[j]);
        double x_diff(obs.x-lm.x_f);
        double y_diff(obs.y-lm.y_f);
        double const distance(sqrt(x_diff*x_diff+y_diff*y_diff));
      if(obs.id==-1)
      {
        obs.id = i;
        obs.distanceTo = distance;
      }
      else
      {
        if(obs.distanceTo>distance)
        {
          obs.id = i;
          obs.distanceTo = distance;
        }
      }
    }
  }
}

void Particle::calculateWeight(double const std_landmark[])
{
  double mean_x(0.), mean_y(0.);
  double const nominator(2*M_PI*std_landmark[0]*std_landmark[1]);
  double const nom2std_x(2*std_landmark[0]*std_landmark[0]);
  double const nom2std_y(2*std_landmark[1]*std_landmark[1]);

  //first step: calculate the mean of x and y

  //MEAN x and MEAN y is the distance of the MapObject and the Observation
  //we finally choose.
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
//  default_random_engine gen;

  double state[3];
  state[0]=x;
  state[1]=y;
  state[2]=theta;

  for(int i(0); i<num_particles; i++)
  {
    enrichWithRandomNoise(state, std);
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
    particles[i].predictPosition(delta_t, std_pos, velocity, yaw_rate);
  }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

  //MScharf: nearest neighboor data association and assign each observation a landmark

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
  for(int i(0); i<num_particles; i++)
  {
    particles[i].transformAndUpdateObservation(observations);
    particles[i].matchObservationsToMap(map_landmarks, sensor_range);
    particles[i].calculateWeight(std_landmarks);
  }

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
//	particle.associations.clear();
//	particle.sense_x.clear();
//	particle.sense_y.clear();
//
//	particle.associations= associations;
// 	particle.sense_x = sense_x;
// 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
//	vector<int> v = best.associations;
//	stringstream ss;
//    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
//    string s = ss.str();
//    s = s.substr(0, s.length()-1);  // get rid of the trailing space
//    return s;
  string s;
  return s;
}
string ParticleFilter::getSenseX(Particle best)
{
//	vector<double> v = best.sense_x;
//	stringstream ss;
//    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
//    string s = ss.str();
//    s = s.substr(0, s.length()-1);  // get rid of the trailing space
//    return s;
  string s;
  return s;
}
string ParticleFilter::getSenseY(Particle best)
{
//	vector<double> v = best.sense_y;
//	stringstream ss;
//    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
//    string s = ss.str();
//    s = s.substr(0, s.length()-1);  // get rid of the trailing space
//    return s;
  string s;
  return s;
}

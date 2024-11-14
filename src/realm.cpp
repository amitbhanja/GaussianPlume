/*
 * uFVM-CPP:
 *      Parallel unstructured 3D fluid flow CFD package 2020
 *
 * Developed by:
 *      Fadl Moukalled
 *      Mhamad Mahdi Alloush
 *      Adam Fares
 *
 * File:
 *      realm.cpp
*/

// Built-in Libraries
#include <chrono>
#include <regex>
#include <thread>
#include <fstream>
#include <sys/stat.h>
#include <unistd.h>
#include <stdexcept>
#include <unordered_map>
#include <iostream>

// Code libraries
#include "realm.h"
//#include "latlon2local.h"

// new libraries
#include "cmath"


namespace ufvm
{


    // Constructors

    realm::realm()
    {
        // Create output directories

    }

    // Operations





    //******* 6.09.2024, Jorian Schlunegger (JS!),  Gaussian Plume  *********************************//

    void realm::run()
    {
        /*
        iter_ = 0;
        timeStepCount_ = 0;
        auto myCentroids = mesh_ref().centroids().internalDataArray();
        
        label numberOfCells = myCentroids.size();
        std::cout << "Number of cells: " << numberOfCells << std::endl << std::endl;

        volField<scalar> C(mesh_ptr(),"Concentration"); 

        //auto C_interior = C.internalDataArray();
        auto U_interior = settings_ref().wind_velocity();
        auto Q = settings_ref().Q_rate();                               // pollution rate in g/s
        //auto u1 = U_interior.magnitude();
        auto passquill_stability = settings_ref().wind_stability();
        Point stack_outlet = settings_ref().stack_outlet_Point();     // Outlet of stack
        Point receiver_point = settings_ref().receiver_point_Point();

        for(label i=0; i < numberOfCells; i++)
        {
            auto x = myCentroids[i].x();
            auto y = myCentroids[i].y();
            auto z = myCentroids[i].z();

            C[i] = compute_Concentration(myCentroids[i], Q, U_interior, passquill_stability, stack_outlet);
        }

        C.write(caseDirectoryPath()/removeTrailingZeros(std::to_string(time_)));

        // ----------------------- This part here: debugging -----------------------//
        // calculate receiver point
        scalar receiver_sigma_y = -1;
        scalar receiver_sigma_z = -1;
        scalar receiver_phi = -1;
        scalar C_receiver_x = 0;
        scalar C_receiver_y = 0;
        scalar receiver_phi_x;

        // calculate wind  
        computeAlpha_x(receiver_phi_x, U_interior);     
        std::cout << "wind speed is " << std::sqrt(pow(U_interior.x(),2) + pow(U_interior.y(),2))<<std::endl;
        std::cout << "x component is " << U_interior.x()<<std::endl;
        std::cout << "y component is " << U_interior.y()<<std::endl;
        std::cout << "angle of wind to x-direction is  " << receiver_phi_x*180/M_PI <<"°"<<std::endl <<std::endl;
                
        Point receiver_oldVect = rotatePointInv(receiver_point, stack_outlet, receiver_phi_x);

        //std::cout << "origin Point receiver would be:" << receiver_oldVect << std::endl;

        if ((receiver_oldVect.x()-stack_outlet.x()>0) )
        {
            computeSigmaYSigmaZ(receiver_oldVect.x()-stack_outlet.x(), passquill_stability, receiver_sigma_y, receiver_sigma_z, receiver_phi);
            
            auto C_receiver = compute_Concentration_old(receiver_oldVect, Q, U_interior, passquill_stability, stack_outlet);

            std::cout<< "Measurement point: "<<receiver_point <<", C = "<< C_receiver<< " g/m³"<<std::endl;
            std::cout<< "sigma y: "<<receiver_sigma_y <<", sigma z: "<<receiver_sigma_z <<", phi: "<<receiver_phi << std::endl << std::endl << std::endl;
        }
        else
        {
            auto C_receiver = 0;
        }
        */
        // ----------------------- end debugging -----------------------//

/*
        // get cartesian coordinates from lat long
        double lat[100];
        double lon[100];
        double origin[3];

        origin[0] = stack_outlet.x();
        origin[1] = stack_outlet.y();
        origin[2] = stack_outlet.z();

        std::cout << "Stack Outlet " << stack_outlet << std::endl;
        coordinates latLongData = receiver_points;

        for(int i = 0; i < receiver_points.size(); ++i)
        {
            auto recv_point = receiver_points[i];
            lat[i] = recv_point.x();
            lon[i] = recv_point.y();
        }
        double x[100];
        double y[100];
        //coder::latlon2local(lat, lon, origin, x, y);

        for(int i = 0; i < receiver_points.size(); ++i)
        {
            std::cout << "Receiver Point " << receiver_points[i] << " (" << x[i] << ", " << y[i] << ")" << std::endl;
            receiver_points[i].x() = x[i];
            receiver_points[i].y() = y[i];
        }
*/

        Point origin_coordinate(0,0,stack_outlet.z());
        for(int i = 0; i < receiver_points.size(); ++i)
        {
            auto recv_point = distance(stack_outlet.x(), stack_outlet.y(), receiver_points[i].x(), receiver_points[i].y());
            std::cout << "Distance of " << receiver_points[i] << " from " << stack_outlet << " " << recv_point << std::endl;
            auto recv_meas = receiver_measurements[i];
            auto wind_vel = wind_velocity[i];
            for(int j = 0; j < 6; ++j)
            {
                scalar conc = compute_Concentration(recv_point, flow_rate, wind_vel, j, origin_coordinate);
                std::string keyword = std::string(receiver_points[i]) + "," + std::to_string(j) + "," + std::to_string(recv_meas);
                forward_model_values.emplace_back(std::make_pair(keyword, std::to_string(conc)));
                std::cout << "Concentration " << recv_point << " with stability class " << j << " : " << conc << std::endl;
            }

            /*
            auto receiver_oldVect = rotatePointInv(recv_point, origin_coordinate, wind_vel.direction());

            for(int j = 0; j < 6; ++j)
            {
                scalar Q_out = compute_amountQ(receiver_oldVect, recv_meas, wind_vel, j, origin_coordinate, {0,0,0}, 0);
                //Q_out *= 3600;
                std::string keyword = std::string(receiver_points[i]) + "," + std::to_string(j);
                reverse_model_values.emplace_back(std::make_pair(keyword, std::to_string(Q_out)));

                std::cout << "Flow Rate for " << recv_point << " with stability class " << j << " : " << Q_out << std::endl;
            }
            */
        }
    }

    int realm::signum(int num) 
    {
        if (num >= 0) 
        {
            return 1;
        } 
        else if (num < 0) 
        {
            return -1;
        } 
        else 
        {
            return 0;
        }
    }

    void realm::computeAlpha_x(scalar& phi_x, Point U)
    {
        if (!((U.x() == 0) && (U.y() == 0)))
        {
            phi_x = signum(U.y())*acos(U.x()/(sqrt(pow(U.x(),2) + pow(U.y(),2))));    // Calculates the angle between Point and x-axis
        }
        else
        {
            phi_x = 0;
        }
        
    }

    Point realm::rotatePointInv(Point& v, Point& stack, scalar angle)
    {   // This is the invert rotation matrix to calculate the values in old coordinate system
        
        // correct difference to stack is taken into account here!
        scalar x = v.x() - stack.x();
        scalar y = v.y() - stack.y();
        //scalar z = v.z() - stack.z();

        scalar cosAngle = cos(angle);
        scalar sinAngle = sin(angle);

        scalar oldX = x * cosAngle + y * sinAngle;
        scalar oldY = -x * sinAngle + y * cosAngle;

        return Point(oldX+stack.x(), oldY+stack.y(), v.z());
    }

    Point realm::distance(scalar lat1, scalar lon1, scalar lat2, scalar lon2)
    {
        constexpr scalar EARTH_RADIUS = 6371000.0;

        double lat1rad = lat1 * M_PI / 180.0;
        double lon1rad = lon1 * M_PI / 180.0;

        double lat2rad = lat2 * M_PI / 180.0;
        double lon2rad = lon2 * M_PI / 180.0;

        scalar x = EARTH_RADIUS * (lon2rad - lon1rad) * cos((lat1rad + lat2rad)/2);
        scalar y = EARTH_RADIUS * (lat2rad - lat1rad);

        return Point(x, y, 0);
    }

    void realm::computeSigmaYSigmaZ(scalar dx, int stability, scalar& sigma_y, scalar& sigma_z, scalar& phi)
    {
        
        scalar a;
        scalar b;
        scalar c;
        scalar d;
        //scalar dx = x-stack_outlet.x();
        //scalar phi;
               
        if (stability == 0)
        {
            c=24.1670;
            d=2.5334;

            if((dx>=0.0) && (dx<100.0))
            {
                a=122.800;
                b=0.94470;
            }

            else if((dx>=100.0) && (dx<150.0))
            {        
                a=158.080;
                b=1.05420;
            }

            else if((dx>=150.0) && (dx<200.0))
            {       
                a=170.220;
                b=1.09320;
            }

            else if((dx>=200.0) && (dx<250.0)) 
            {       
                a=179.520;
                b=1.12620;
            }

            else if((dx>=250.0) && (dx<300.0)) 
            {       
                a=217.410;
                b=1.26440;
            }

            else if((dx>=300.0) && (dx<400.0))
            {        
                a=258.890;
                b=1.40940;
            }

            else if((dx>=400.0) && (dx<500.0))  
            {      
                a=346.750;
                b=1.72830;
            }

            else if((dx>=500.0) && (dx<3110.0))
            {        
                a=453.850;
                b=2.11660;
            }

            else if(dx>=3110.0)
            {
                a=453.850;
                b=2.11660;
            }

            else
            {
                a=1.0;
                b=1.0;
            }
        }
    
        else if (stability == 1)
        {
            c=18.3330;
            d=1.8096;

            if((dx>=0.0) && (dx<200.0))
            {
                a=90.673;
                b=0.93198;
            }

            else if((dx>=200.0) && (dx<400.0))
            {        
                a=98.483;
                b=0.98332;
            }

            else if(dx>=400.0)
            {
                a=109.300;
                b=1.09710;
            }

            else
            {
                a=1.0;
                b=1.0;
            }
        }

        else if (stability == 2)
        {
            a=61.141;
            b=0.91465;
            c=12.5000;
            d=1.0857;
        }
            
        else if (stability == 3)
        {
            c=8.3330;
            d=0.72382;

            if ((dx>=0.0) && (dx<300))
            {
                a=34.459;
                b=0.86974;
            }

            else if ((dx>=300.0) && (dx<1000.0))
            {
                a=32.093;
                b=0.81066;
            }

            else if ((dx>=1000.0) && (dx<3000.0))
            {
                a=32.093;
                b=0.64403;
            }

            else if ((dx>=1000.0) && (dx<10000.0))
            {
                a=33.504;
                b=0.60486;
            }

            else if ((dx>=10000.0) && (dx<30000.0))
            {
                a=36.650;
                b=0.56589;
            }

            else if (dx>=30000.0)
            {
                a=44.053;
                b=0.51179;
            }

            else
            {
                a=1;
                b=1;
            }
        }

        else if (stability == 4)
        {
            c=6.2500;
            d=0.54287;

            if ((dx>=0.0) && (dx<100))
            {
                a=24.260;
                b=0.83660;
            }

            else if ((dx>=100.0) && (dx<300.0))
            {
                a=23.331;
                b=0.81956;
            }

            else if ((dx>=300.0) && (dx<1000.0))
            {
                a=21.628;
                b=0.75660;
            }

            else if ((dx>=1000.0) && (dx<2000.0))
            {
                a=21.628;
                b=0.63077;
            }

            else if ((dx>=2000.0) && (dx<4000.0))
            {
                a=22.534;
                b=0.57154;
            }

            else if ((dx>=4000.0) && (dx<10000.0))
            {
                a=24.703;
                b=0.50527;
            }

            else if ((dx>=10000.0) && (dx<20000.0))
            {
                a=26.970;
                b=0.46713;
            }

            else if ((dx>=20000.0) && (dx<40000.0))
            {
                a=35.420;
                b=0.37615;
            }

            else if (dx>=40000.0)
            {
                a=47.618;
                b=0.29592;
            }

            else
            {
                a=1;
                b=1;
            }
        }

        else if (stability == 5)
        {
            c=4.1667;
            d=0.36191;

            if ((dx>0.0) && (dx<200))
            {
                a=15.209;
                b=0.81558;
            }

            else if ((dx>200) && (dx<700))
            {
                a=14.457;
                b=0.78407;
            }

            else if ((dx>700) && (dx<1000))
            {
                a=13.953;
                b=0.68465;
            }

            else if ((dx>1000) && (dx<2000))
            {
                a=13.953;
                b=0.63227;
            }

            else if ((dx>2000) && (dx<3000))
            {
                a=14.823;
                b=0.54503;
            }

            else if ((dx>3000) && (dx<7000))
            {
                a=16.187;
                b=0.46490;
            }

            else if ((dx>7000) && (dx<15000))
            {
                a=17.836;
                b=0.41507;
            }

            else if ((dx>15000) && (dx<30000))
            {
                a=22.651;
                b=0.32681;
            }

            else if ((dx>30000) && (dx<60000))
            {
                a=27.074;
                b=0.27436;
            }

            else if (dx>60000)
            {
                a=34.219;
                b=0.21716;
            }

            else
            {
                a=1;
                b=1;
            }
        }

        else
        {
            a=1;
            b=1;
            c=1;
            d=1;
        }

        // dx has to be in (km) for following equations
        sigma_z = a*pow(dx/1000,b);     
        if (sigma_z > 5000)
        {
            sigma_z = 5000;
        }
            
        phi = 0.017453293*(c-d*log(dx/1000));
        sigma_y = 465.11628*dx/1000*tan(phi);
    }

    scalar realm::compute_Concentration(Point cell, scalar Q, Velocity U, int stabilityClass, Point& stack_outlet)
    {   
        scalar windAngle = U.direction();
        
        Point oldVect = rotatePointInv(cell, stack_outlet, windAngle);

        scalar x = oldVect.x();
        scalar y = oldVect.y();
        scalar z = cell.z();

        scalar U_mag = U.magnitude();
        //scalar U_mag = std::sqrt(pow(U.x(),2) + pow(U.y(),2));

        scalar sigma_y = -1;
        scalar sigma_z = -1;
        scalar phi = -1;
        scalar dx = x-stack_outlet.x();
        
        if (dx>0)  // Concentration value is only calculated for cell after the stack in x direction
        {
            computeSigmaYSigmaZ(dx, stabilityClass, sigma_y, sigma_z, phi);
            return Q /(2.0 * M_PI * U_mag * sigma_y * sigma_z) * exp(-pow((y-stack_outlet.y()), 2) / (2.0 * pow(sigma_y, 2))) * 
                (exp(-pow(z - stack_outlet.z(), 2) / (2.0 * pow(sigma_z, 2))) + exp(-pow(z + stack_outlet.z(), 2) / (2.0 * pow(sigma_z, 2))));
        }
        
        else
        {
            return 0;
        }
    }

    scalar realm::compute_amountQ(Point& cell, scalar C_meas, Velocity& U, int stabilityClass, Point& stack_outlet, Point GPS_unc, scalar Sensor_unc)
    {
        scalar sigma_y = -1;
        scalar sigma_z = -1;
        scalar phi = -1;
        scalar uncX = 1 + GPS_unc.x();
        scalar uncY = 1 + GPS_unc.y();
        scalar uncZ = 1 + GPS_unc.z();
        scalar uncSens = 1 + Sensor_unc;
        std::cout << "Cell " << cell << std::endl;
        std::cout << "Stack Outlet " << stack_outlet << std::endl;

        std::cout << "Diff " << cell.x() - stack_outlet.x() << std::endl;

        if ((cell.x() - stack_outlet.x())>0)  // Concentration value is only calculated for cell after the stack in x direction
        {
            computeSigmaYSigmaZ((cell.x() - stack_outlet.x())*uncX, stabilityClass, sigma_y, sigma_z, phi);

            return C_meas*uncSens * (2.0 * M_PI * U.magnitude() * sigma_y * sigma_z) / ( exp(-pow((cell.y()-stack_outlet.y())*uncY, 2) / (2.0 * pow(sigma_y, 2))) * 
                (exp(-pow((cell.z() - stack_outlet.z())*uncZ, 2) / (2.0 * pow(sigma_z, 2))) + exp(-pow((cell.z() + stack_outlet.z())*uncZ, 2) / (2.0 * pow(sigma_z, 2)))));
          
        }
        else
        {
            return 0;
        }
    }

    // Following function is used for debugging
    scalar realm::compute_Concentration_old(Point cell, scalar Q, Point U, int& stabilityClass, Point& stack_outlet)
    {   
        scalar x = cell.x();
        scalar y = cell.y();
        scalar z = cell.z();

        scalar sigma_y = -1;
        scalar sigma_z = -1;
        scalar phi = -1;
        scalar dx = x-stack_outlet.x();
        //scalar phi_x = atan(U.y()/(U.x() + 1e-12));
        //scalar dx = cos(phi_x)*(x-stack_outlet.x());
        //scalar dy = sin(phi_x)*(y-stack_outlet.y());
        
        if (dx>=0)  // Concentration value is only calculated for cell after the stack in x direction
        {
            computeSigmaYSigmaZ(dx, stabilityClass, sigma_y, sigma_z, phi);
            return Q /(2.0 * M_PI * U.magnitude() * sigma_y * sigma_z) * exp(-pow((y-stack_outlet.y()), 2) / (2.0 * pow(sigma_y, 2))) * 
                (exp(-pow(z - stack_outlet.z(), 2) / (2.0 * pow(sigma_z, 2))) + exp(-pow(z + stack_outlet.z(), 2) / (2.0 * pow(sigma_z, 2))));
        }

        else
        {
            return 0;
        }
    }
    
    // ********************************************************************************************//

} // accel

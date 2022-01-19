/*

Copyright (c) 2005-2022, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifdef CHASTE_CVODE
#ifndef _ABSTRACTCVODECELLWITHDATACLAMP_HPP_
#define _ABSTRACTCVODECELLWITHDATACLAMP_HPP_

#include "AbstractCvodeCell.hpp"

/**
 * A cardiac cell that is designed to be simulated using CVODE.
 * It uses CVODE's vector type natively via AbstractCvodeSystem.
 *
 * This subclass of AbstractCvodeCell contains methods used by systems that have
 * had an 'artificial' data clamp current introduced into the dV/dt equation.
 * These are generated by chaste_codegen with ConvertCellModel using the option '--cvode-data-clamp'.
 *
 */
class AbstractCvodeCellWithDataClamp : public AbstractCvodeCell
{
private:
    /** For testing the interpolation method */
    friend class TestCvodeCellsWithDataClamp;

    /** Needed for serialization. */
    friend class boost::serialization::access;

    /**
     * Archive the member variables.
     *
     * @param archive
     * @param version
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // This calls serialize on the base class.
        archive & boost::serialization::base_object<AbstractCvodeCell>(*this);

        // Archive this class' member variables.
        archive & mExperimentalTimes;
        archive & mExperimentalVoltages;
        archive & mDataClampIsOn;
        archive & mDataAvailable;
    }

    /** The experimental times of voltage measurements in the data */
    std::vector<double> mExperimentalTimes;

    /** The voltages recorded in experiment at the times in #mExperimentalTimes, forming data for clamp. */
    std::vector<double> mExperimentalVoltages;

    /** Whether experimental data have been set */
    bool mDataAvailable;

protected:

    /** Whether the data clamp is active at the moment */
    bool mDataClampIsOn;

    /**
     * Linear interpolation method to work out the voltage at any given time the ODE solver requests.
     *
     * @param rTime  The time in the data trace at which to return the voltage (linearly interpolated).
     * @return The voltage at rTime.
     */
    double GetExperimentalVoltageAtTimeT(const double& rTime);

public:
    /**
     * Constructor
     *
     * @param pSolver  UNUSED for CVODE cells.
     * @param numberOfStateVariables  Number of ODEs in the system.
     * @param voltageIndex  The index of the variable which represents membrane voltage.
     * @param pIntracellularStimulus  The intracellular stimulus function.
     */
    AbstractCvodeCellWithDataClamp(boost::shared_ptr<AbstractIvpOdeSolver> pSolver,
                                   unsigned numberOfStateVariables,
                                   unsigned voltageIndex,
                                   boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus);

    /**
     * Destructor: free any memory we've used.
     */
    virtual ~AbstractCvodeCellWithDataClamp();

    /**
     * Provide experimental data. Must be called before #TurnOnDataClamp().
     *
     * @param experimentalTimes  a vector of times - must be monotonically increasing, but can have variable intervals.
     * @param experimentalVoltages  a vector of voltages recorded at these times.
     *
     * These two vectors must be the same length.
     */
    void SetExperimentalData(std::vector<double> experimentalTimes,
                             std::vector<double> experimentalVoltages);


    /**
     * Switch off the data clamping current
     */
    void TurnOffDataClamp();

    /**
     * Switch on the data clamping current
     *
     * @param conductance  The conductance of the data clamping current - don't yet know best number to set this to.
     */
    void TurnOnDataClamp(double conductance = 100);
};

CLASS_IS_ABSTRACT(AbstractCvodeCellWithDataClamp)

#endif // _ABSTRACTCVODECELLWITHDATACLAMP_HPP_
#endif // CHASTE_CVODE

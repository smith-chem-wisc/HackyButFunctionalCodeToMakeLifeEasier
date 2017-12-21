using System;
using System.Collections.Generic;

namespace FusionPeptideMassVerification
{
    public class MS2Spectra
    {
        public static readonly int decimalsToRound = 5;
        public readonly int scanNumber;
        public readonly double precursorMass;
        public readonly double[] masses;
        public readonly double[] intensities;
        public readonly double totalIntensity;

        public MS2Spectra(int scanNumber, double precursorMass, double[] mzs, double[] intensities)
        {
            this.scanNumber = scanNumber;
            this.precursorMass = precursorMass;
            List<double> masses = new List<double> { 0 };
            List<double> intensitiesList = new List<double> { 1 };
            for (int i = 0; i < mzs.Length; i++)
            {
                masses.Add(Math.Round(mzs[i] - Constants.PROTON_MASS, decimalsToRound));
                double roundedIntensity = Math.Round(intensities[i], decimalsToRound);
                intensitiesList.Add(roundedIntensity);
                totalIntensity += roundedIntensity;
            }
            masses.Add(precursorMass);
            intensitiesList.Add(1);
            this.intensities = intensitiesList.ToArray();
            this.masses = masses.ToArray();
        }
    }
}

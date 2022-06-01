using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using System.Text;
using System.Threading.Tasks;
using IsotopomerDynamics;
using MathNet.Numerics.Statistics;

namespace RsquaredAndPearsonComputationForMassIsotopmers
{
    public class RsquaredComputation
    {
        public static double ph = 1.5574E-4;
        public RsquaredComputation()
        {

        }

        public List<double> computation(float[] fTheoreticalIsotopes, List<float> fTimeCourse, float szNEH,
            Object fBodyWaterEnrichment, Object TimeCourseIsotopeCluster,
            Object aTimePointsForThisPeptde, string szRateFile, float Rateconst, float new_RMSE)
        {
            List<List<float>> experimentslist = (List<List<float>>)TimeCourseIsotopeCluster;
            List<ExperimentRecord> experiments = new List<ExperimentRecord>();
            var temp_fBodyWaterEnrichment = ((List<float>)fBodyWaterEnrichment);
            float M0 = 100 * fTheoreticalIsotopes[0] / (fTheoreticalIsotopes[0] + fTheoreticalIsotopes[1] +
                fTheoreticalIsotopes[2] + fTheoreticalIsotopes[3] + fTheoreticalIsotopes[4] + fTheoreticalIsotopes[5]);

            float final_BWE = temp_fBodyWaterEnrichment[temp_fBodyWaterEnrichment.Count - 1];


            for (int i = 0; i < experimentslist.Count; i++)
            {
                ExperimentRecord er = new ExperimentRecord();
                er.I0 = experimentslist[i][0];
                er.I1 = experimentslist[i][1];
                er.I2 = experimentslist[i][2];
                er.I3 = experimentslist[i][3];
                er.I4 = experimentslist[i][4];
                er.I5 = experimentslist[i][5];
                er.ExperimentTime = fTimeCourse[i];// (double)fScaleTime;
                er.Deuteriumenrichment = temp_fBodyWaterEnrichment[i];
                experiments.Add(er);
            }


            experiments = computeDeuteriumenrichmentInPeptide((double)szNEH, fTheoreticalIsotopes, szRateFile, experiments, M0);
            var RIAvalues = computeRIAPerExperiment(experiments);
            var normalizedValues = normalizeRIAValuesForAllPeptides(RIAvalues, M0, final_BWE, szNEH, temp_fBodyWaterEnrichment);

            List<RIA> mergedRIAvalues = mergeMultipleRIAPerDay2(normalizedValues, fTimeCourse);

            List<TheoreticalI0Value> theoreticalI0Values = computeTheoreticalCurvePoints(fTimeCourse, final_BWE, Rateconst, M0, szNEH);

            double newPearsonCorrval = computreNewPearsonCorrlationValue(mergedRIAvalues.Select(x => x.I0_t_fromA1A0).ToList(),
                             mergedRIAvalues.Select(x => x.I0_t_fromA2A0).ToList(),
                             mergedRIAvalues.Select(x => x.I0_t_fromA2A1).ToList(),
                             mergedRIAvalues.Select(x => x.RIA_value).ToList(),
                            theoreticalI0Values.Select(x => x.value).ToList(), ref new_RMSE);

            return new List<double>() { newPearsonCorrval, new_RMSE };

        }

        public double computreNewPearsonCorrlationValue(List<double?> a1ao, List<double?> a2ao, List<double?> a1a2,
            List<double?> experimental_RIA, List<double> theoretical_RIA, ref float new_RMSE)
        {
            try
            {
                var selected_points = new List<double>();
                var selected_A1A0_count = 0;
                var selected_A2A0_count = 0;
                var selected_A2A1_count = 0;
                for (int i = 0; i < experimental_RIA.Count; i++)
                {
                    var theoretical_val = (double)theoretical_RIA[i];
                    var candidate_points = new List<double>();

                    candidate_points.Add(a1ao[i] == null || double.IsNaN((double)a1ao[i]) ? double.MaxValue : Math.Abs((double)a1ao[i] - theoretical_val));
                    candidate_points.Add(a2ao[i] == null || double.IsNaN((double)a2ao[i]) ? double.MaxValue : Math.Abs((double)a2ao[i] - theoretical_val));
                    candidate_points.Add(a1a2[i] == null || double.IsNaN((double)a1a2[i]) ? double.MaxValue : Math.Abs((double)a1a2[i] - theoretical_val));
                    candidate_points.Add(experimental_RIA[i] == null || double.IsNaN((double)experimental_RIA[i]) ? double.MaxValue : Math.Abs((double)experimental_RIA[i] - theoretical_val));

                    // index of minimum point
                    var min_val = candidate_points.Min();

                    // add the minimum error point to the selected list for the specific time point
                    if (min_val == double.MaxValue) selected_points.Add(double.NaN);
                    else
                    {
                        var index_min_val = candidate_points.IndexOf(min_val);
                        switch (index_min_val)
                        {
                            case 0: selected_points.Add((double)a1ao[i]); selected_A1A0_count += 1; break;
                            case 1: selected_points.Add((double)a2ao[i]); selected_A2A0_count += 1; break;
                            case 2: selected_points.Add((double)a1a2[i]); selected_A2A1_count += 1; break;
                            case 3: selected_points.Add((double)experimental_RIA[i]); break;
                            default: selected_points.Add(double.NaN); break;
                        }
                    }
                }



                var rsquared = computeRsquared(selected_points, theoretical_RIA);
                var rsquaredold = computeRsquared(experimental_RIA.Select(x => (double)x).ToList(), theoretical_RIA);

                // filter null values from experimental values
                List<List<double>> temp = new List<List<double>>();
                temp.Add(new List<double>()); // experimental value
                temp.Add(new List<double>()); // theoretical value
                temp.Add(new List<double>()); // selected value

                for (int i = 0; i < experimental_RIA.Count; i++)
                {
                    if (!double.IsNaN(experimental_RIA[i].Value))
                    {
                        temp[0].Add(experimental_RIA[i].Value);
                        temp[1].Add(theoretical_RIA[i]);
                        temp[2].Add(selected_points[i]);
                    }
                }

                var _correlation_old = Correlation.Pearson(temp[0], temp[1]);
                var _correlation = Correlation.Pearson(temp[2], temp[1]);

                var correlation_old = PearsonCorrelation(temp[0], temp[1]);
                var correlation = PearsonCorrelation(temp[2], temp[1]);

                new_RMSE = (float)computeRMSE(temp[2], temp[1]);


                Console.WriteLine("===>");
                Console.WriteLine(String.Format("Old Pearson Correlation = {0}, New Pearson Correlation = {1}, New RMSE = {2}",
                    correlation_old, correlation, new_RMSE));

                return correlation;

            }
            catch (Exception ex)
            {
                Console.WriteLine(ex.ToString());
            }

            return double.NaN;
        }

        public double PearsonCorrelation(List<double> FirstArr, List<double> SecondArr)
        {
            double fMean1, fMean2, ftemp, fsquare1, fsquare2;
            int nPoints = FirstArr.Count;
            int k;



            fMean1 = fMean2 = ftemp = 0;

            for (k = 0; k < nPoints; k++)
            {
                fMean1 = fMean1 + FirstArr[k];

                fMean2 = fMean2 + SecondArr[k];
            }

            fMean1 = fMean1 / (float)nPoints;

            fMean2 = fMean2 / (float)nPoints;

            fsquare1 = fsquare2 = 0;

            for (k = 0; k < nPoints; k++)
            {
                fsquare1 = fsquare1 + (fMean1 - FirstArr[k]) * (fMean1 - FirstArr[k]);

                fsquare2 = fsquare2 + (fMean2 - SecondArr[k]) * (fMean2 - SecondArr[k]);

                ftemp = ftemp + (fMean1 - FirstArr[k]) * (fMean2 - SecondArr[k]);
            }

            if (fsquare1 < 0.00000001 || fsquare2 < 0.00000001)
                return -2.0f;

            ftemp = ftemp / Math.Sqrt(fsquare1 * fsquare2);

            return ftemp;
        }

        public double computeRsquared(List<double> experimentalValue, List<double> fitvalue)
        {
            var mean_exp = experimentalValue.Where(x => !double.IsNaN(x)).Average();
            double rss = 0;
            double ss = 0;
            double rsquared = double.NaN;

            for (int i = 0; i < experimentalValue.Count; i++)
            {
                if (!double.IsNaN(experimentalValue[i]))
                {
                    ss = ss + Math.Pow((double)(experimentalValue[i] - mean_exp), 2);
                    rss = rss + Math.Pow((double)(experimentalValue[i] - fitvalue[i]), 2);
                }
            }

            //if (r.Rateconst > 0.0006) RSquare = 1 - (rss / ss);
            //else RSquare = 1 - (diff);

            if (ss != 0)
                rsquared = 1 - (rss / ss);

            return rsquared;
        }

        public double computeRMSE(List<double> experimentalValue, List<double> fitvalue)
        {
            var mean_exp = experimentalValue.Where(x => !double.IsNaN(x)).Average();
            double rss = 0;
            double ss = 0;
            double rsquared = double.NaN;

            for (int i = 0; i < experimentalValue.Count; i++)
            {
                if (!double.IsNaN(experimentalValue[i]))
                {
                    ss = ss + Math.Pow((double)(experimentalValue[i] - mean_exp), 2);
                    rss = rss + Math.Pow((double)(experimentalValue[i] - fitvalue[i]), 2);
                }
            }

            var RMSE_value = Math.Sqrt(rss / experimentalValue.Count());

            return RMSE_value;
        }

        public List<TheoreticalI0Value> computeTheoreticalCurvePoints(List<float> fTimeCourse, float final_BWE,
            float Rateconst, float M0, float NEH)
        {
            List<TheoreticalI0Value> theoreticalI0Values = new List<TheoreticalI0Value>();
            try
            {
                double ph = 1.5574E-4;
                double pw = final_BWE;
                double io = 0;
                double neh = 0;
                double k = 0;

                {


                    List<double> mytimelist = new List<double>();
                    var experiment_time = fTimeCourse.Select(x => x).Distinct().ToList();

                    foreach (double t in experiment_time)
                    {
                        //foreach (int t in this.Experiment_time)
                        //{
                        try
                        {
                            io = (double)(M0 / 100);
                            neh = (double)(NEH);
                            k = (double)(Rateconst);

                            var val1 = io * Math.Pow(1 - (pw / (1 - ph)), neh);
                            var val2 = io * Math.Pow(Math.E, -1 * k * t) * (1 - (Math.Pow(1 - (pw / (1 - ph)), neh)));


                            var val = val1 + val2;

                            TheoreticalI0Value theoreticalI0Value = new TheoreticalI0Value("", t, val, 0);
                            theoreticalI0Values.Add(theoreticalI0Value);
                        }
                        catch (Exception ex)
                        {
                            Console.WriteLine("Error => computeExpectedCurvePoints(), " + ex.Message);
                            continue;
                        }

                    }
                }
            }
            catch (Exception e) { Console.WriteLine("Error => computeExpectedCurvePoints(), " + e.Message); }
            return theoreticalI0Values;
        }

        public List<RIA> mergeMultipleRIAPerDay2(List<RIA> RIAvalues, List<float> fTimeCourse)
        {
            List<RIA> temp_RIAvalues = RIAvalues;
            List<RIA> mergedRIAvalues = new List<RIA>();
            var distincttime = fTimeCourse.Select(x => x).Distinct().ToList();
            foreach (float t in distincttime)
            {
                List<RIA> temp_RIAvalues_pertime = temp_RIAvalues.Where(x => Math.Abs(x.Time - t) <= 0.001).ToList();

                // check the number of experiments with zero ion score
                int countOfNonZeroIonScore = temp_RIAvalues_pertime.Where(x => x.IonScore > 0).Count();
                if (countOfNonZeroIonScore > 0) temp_RIAvalues_pertime = temp_RIAvalues_pertime.Where(x => x.IonScore > 0).ToList();

                RIA ria = new RIA();
                ria.ExperimentNames = new List<string>();
                ria.Time = t;

                double sum_io = 0; for (int i = 0; i < temp_RIAvalues_pertime.Count(); i++) sum_io = (double)(sum_io + temp_RIAvalues_pertime[i].I0);
                double sum_ioria = 0; for (int i = 0; i < temp_RIAvalues_pertime.Count(); i++) sum_ioria = (double)(
                        sum_ioria + (temp_RIAvalues_pertime[i].I0 * temp_RIAvalues_pertime[i].RIA_value));
                var new_ria = sum_ioria / sum_io;

                #region compute the modified I0_t

                ria.I0_t_fromA1A0 = temp_RIAvalues_pertime.Count > 0 ? temp_RIAvalues_pertime.Select(x => x.I0_t_fromA1A0).Average() : null;
                ria.I0_t_fromA2A0 = temp_RIAvalues_pertime.Count > 0 ? temp_RIAvalues_pertime.Select(x => x.I0_t_fromA2A0).Average() : null;
                ria.I0_t_fromA2A1 = temp_RIAvalues_pertime.Count > 0 ? temp_RIAvalues_pertime.Select(x => x.I0_t_fromA2A1).Average() : null;
                #endregion

                ria.RIA_value = new_ria;
                ria.ExperimentNames = temp_RIAvalues_pertime.Select(x => x.ExperimentName).ToList();
                mergedRIAvalues.Add(ria);
            }

            return mergedRIAvalues;
        }

        public List<RIA> normalizeRIAValuesForAllPeptides(List<RIA> datapoints, float M0, float final_bwe, float neh,
            List<float> fBodyWaterEnrichment)
        {
            List<RIA> normalizedRIAvalues = new List<RIA>();

            var I0 = 0.0;
            var zeroIimePoints = datapoints.Where(x => x.Time == 0 && x.RIA_value != null && x.RIA_value > 0).ToList();
            if (zeroIimePoints.Count > 0)
                I0 = (double)zeroIimePoints.Select(x => x.RIA_value).Average();

            var temp_mo = (double)M0 / 100;
            if (Math.Abs(I0 - temp_mo) / temp_mo > 0.1) { I0 = temp_mo; }

            var IO_asymptote = I0 * Math.Pow(1 - (final_bwe / (1 - ph)), neh);

            //foreach (var datapoint in datapoints)
            for (int i = 0; i < datapoints.Count; i++)

            {
                var datapoint = datapoints[i];
                //var BWE_t = filecontents.Where(x => x.experimentID == datapoint.ExperimentName).Select(x => x.BWE).FirstOrDefault();
                var BWE_t = fBodyWaterEnrichment[i];


                if (datapoint.Time != 0 && BWE_t == 0) continue;


                var IO_t_asymptote = I0 * Math.Pow(1 - (BWE_t / (1 - ph)), neh);

                var I0_t = (datapoint.Time != 0) ? IO_asymptote + (datapoint.RIA_value - IO_t_asymptote) / (I0 - IO_t_asymptote) * (I0 - IO_asymptote) : datapoint.RIA_value;

                datapoint.RIA_value = I0_t;

                datapoint.I0_t_fromA1A0 = (datapoint.Time != 0) ? IO_asymptote + (datapoint.I0_t_fromA1A0 - IO_t_asymptote) / (I0 - IO_t_asymptote) * (I0 - IO_asymptote) : datapoint.I0_t_fromA1A0;
                datapoint.I0_t_fromA2A0 = (datapoint.Time != 0) ? IO_asymptote + (datapoint.I0_t_fromA2A0 - IO_t_asymptote) / (I0 - IO_t_asymptote) * (I0 - IO_asymptote) : datapoint.I0_t_fromA2A0;
                datapoint.I0_t_fromA2A1 = (datapoint.Time != 0) ? IO_asymptote + (datapoint.I0_t_fromA2A1 - IO_t_asymptote) / (I0 - IO_t_asymptote) * (I0 - IO_asymptote) : datapoint.I0_t_fromA2A1;

                normalizedRIAvalues.Add(datapoint);
            }


            return normalizedRIAvalues;

        }

        public List<ExperimentRecord> computeDeuteriumenrichmentInPeptide(double NEH, Object fTheoreticalIsotopess, string PeptideSeq,
            List<ExperimentRecord> experimentRecordsPerPeptide, float M0)
        {
            var fTheoreticalIsotopes = ((Single[])fTheoreticalIsotopess);


            if (experimentRecordsPerPeptide.Count > 0)
            {
                float[] fNatIsotopes = new float[10];
                float[] fLabIsotopes = new float[10];

                fNatIsotopes[0] = (float)((double)fTheoreticalIsotopes[0]);
                fNatIsotopes[1] = (float)((double)fTheoreticalIsotopes[1]);
                fNatIsotopes[2] = (float)((double)fTheoreticalIsotopes[2]);
                fNatIsotopes[3] = (float)((double)fTheoreticalIsotopes[3]);
                fNatIsotopes[4] = (float)((double)fTheoreticalIsotopes[4]);
                fNatIsotopes[5] = (float)((double)fTheoreticalIsotopes[5]);

                MassIsotopomers MIDyn = new MassIsotopomers();
                int Nall_Hyd = MIDyn.NumberOfHydrogens(PeptideSeq);

                List<float> theo_a2 = new List<float>();
                List<float> theo_a1 = new List<float>();
                List<float> theo_a0 = new List<float>();
                List<double> pxts = new List<double>();

                // compute theoretical values
                for (double bwe = 0.0001; bwe < 0.05; bwe = bwe + 0.0005)
                {
                    MIDyn.CalculateMIDynamics(fNatIsotopes, fLabIsotopes, (float)bwe, (float)NEH, Nall_Hyd);

                    var new_a2_t = fLabIsotopes[2];
                    var new_a1_t = fLabIsotopes[1];
                    var new_a0_t = fLabIsotopes[0];

                    theo_a2.Add(new_a2_t);
                    theo_a1.Add(new_a1_t);
                    theo_a0.Add(new_a0_t);
                    pxts.Add(bwe);
                }


                foreach (ExperimentRecord er in experimentRecordsPerPeptide)
                {
                    double experiment_BWE = (double)er.Deuteriumenrichment;

                    var experimentsAt_t = experimentRecordsPerPeptide.Where(t => t.ExperimentTime == er.ExperimentTime && t.I0 != null && t.I0 > 0 &&
                    (er.I0 + er.I1 + er.I2 + er.I3 + er.I4 + er.I5) > 0
                    ).ToList();
                    if (experimentsAt_t.Count == 0) continue;

                    double sum_A2_t = experimentsAt_t.Sum(x => x.I2).Value;
                    double sum_A2_r_t = experimentsAt_t.Sum(x => (x.I2 * (x.I2 / (x.I0 + x.I1 + x.I2 + x.I3 + x.I4 + x.I5)))).Value;
                    double exp_A2 = sum_A2_r_t / sum_A2_t;

                    double sum_A1_t = experimentsAt_t.Sum(x => x.I1).Value;
                    double sum_A1_r_t = experimentsAt_t.Sum(x => (x.I1 * (x.I1 / (x.I0 + x.I1 + x.I2 + x.I3 + x.I4 + x.I5)))).Value;
                    double exp_A1 = sum_A1_r_t / sum_A1_t;

                    double sum_A0_t = experimentsAt_t.Sum(x => x.I0).Value;
                    double sum_A0_r_t = experimentsAt_t.Sum(x => (x.I0 * (x.I0 / (x.I0 + x.I1 + x.I2 + x.I3 + x.I4 + x.I5)))).Value;
                    double exp_A0 = sum_A0_r_t / sum_A0_t;

                    // compute al/a0
                    if (exp_A0 > 0)
                    {
                        var exp_a1_a0 = exp_A1 / exp_A0;
                        var diff = theo_a1.Select((value, index) => Math.Abs((value / theo_a0[index]) - exp_a1_a0)).ToList();
                        double min_value = diff.Min();
                        int index_min_value = diff.IndexOf(min_value);
                        var selected_pxt = pxts[index_min_value];
                        er.I0_t_fromA1A0 = (double)((M0 / 100.0) * Math.Pow((double)(1 - (selected_pxt / (1 - ph))), (double)NEH));
                        //er.I0_t_fromA1A0_pxt = experiment_BWE;
                    }

                    // compute a2/a0
                    if (exp_A0 > 0)
                    {
                        var exp_a2_a0 = exp_A2 / exp_A0;
                        var diff = theo_a2.Select((value, index) => Math.Abs((value / theo_a0[index]) - exp_a2_a0)).ToList();
                        double min_value = diff.Min();
                        int index_min_value = diff.IndexOf(min_value);
                        var selected_pxt = pxts[index_min_value];
                        er.I0_t_fromA2A0 = (double)((M0 / 100.0) * Math.Pow((double)(1 - (selected_pxt / (1 - ph))), (double)NEH));
                        //er.I0_t_fromA2A0_pxt = experiment_BWE;
                    }


                    // compute a2/a1
                    if (exp_A0 > 0)
                    {
                        var exp_a2_a1 = exp_A2 / exp_A1;
                        var diff = theo_a2.Select((value, index) => Math.Abs((value / theo_a1[index]) - exp_a2_a1)).ToList();
                        double min_value = diff.Min();
                        int index_min_value = diff.IndexOf(min_value);
                        var selected_pxt = pxts[index_min_value];
                        er.I0_t_fromA2A1 = (double)((M0 / 100.0) * Math.Pow((double)(1 - (selected_pxt / (1 - ph))), (double)NEH));
                        //er.I0_t_fromA2A1_pxt = experiment_BWE;
                    }


                }

            }

            return experimentRecordsPerPeptide;
        }

        public List<RIA> computeRIAPerExperiment(List<ExperimentRecord> experimentRecords)
        {
            List<RIA> RIAvalues = new List<RIA>();

            foreach (ExperimentRecord er in experimentRecords)
            {

                var tempsum = er.I0 + er.I1 + er.I2 + er.I3 + er.I4 + er.I5;
                double sum_val = tempsum != null ? (double)(er.I0 + er.I1 + er.I2 + er.I3 + er.I4 + er.I5) : 0;
                if (sum_val != 0)
                {
                    RIA ria = new RIA();
                    ria.ExperimentName = er.ExperimentName;
                    ria.Time = er.ExperimentTime;
                    ria.PeptideSeq = er.PeptideSeq;
                    ria.RIA_value = er.I0 / sum_val;

                    ria.I0 = er.I0;
                    ria.Charge = er.Charge;
                    ria.IonScore = er.IonScore;

                    ria.I0_t_fromA1A0 = er.I0_t_fromA1A0;
                    ria.I0_t_fromA2A0 = er.I0_t_fromA2A0;
                    ria.I0_t_fromA2A1 = er.I0_t_fromA2A1;

                    RIAvalues.Add(ria);
                }
            }
            return RIAvalues;
        }

        public class ExperimentRecord
        {
            public string ExperimentName { get; set; }
            public double ExperimentTime { get; set; }
            public string PeptideSeq { get; set; }
            public double? Charge { get; set; }
            public double SpecMass { get; set; }
            public double? IonScore { get; set; }
            public double? Expectn { get; set; }
            public double? Error { get; set; }
            public double? Scan { get; set; }
            public double? I0 { get; set; }
            public double? I1 { get; set; }
            public double? I2 { get; set; }
            public double? I3 { get; set; }
            public double? I4 { get; set; }
            public double? I5 { get; set; }
            public double? Start_Elution { get; set; }
            public double? End_Elution { get; set; }
            public double? I0_Peak_Width { get; set; }
            public double? Total_Labeling { get; set; }
            public double? Net_Labeling { get; set; }

            public double? Deuteriumenrichment { get; set; } //pX(t) 
            public double? I0_t_fromA1A0 { get; set; }
            public double? I0_t_fromA2A0 { get; set; }
            public double? I0_t_fromA2A1 { get; set; }

            //public double? I0_t_fromA1A0_pxt { get; set; }
            //public double? I0_t_fromA2A0_pxt { get; set; }
            //public double? I0_t_fromA2A1_pxt { get; set; }

        }
        public class RIA
        {
            public string ExperimentName { get; set; }
            public List<string> ExperimentNames { get; set; }
            public string PeptideSeq { get; set; }
            public double? Charge { get; set; }
            public double Time { get; set; }
            public double? I0 { get; set; }
            public double? IonScore { get; set; }
            public double? RIA_value { get; set; }
            public double? I0_t_fromA1A0 { get; set; }
            public double? I0_t_fromA2A0 { get; set; }
            public double? I0_t_fromA2A1 { get; set; }
            //public double? I0_t_fromA1A0_pxt { get; set; }
            //public double? I0_t_fromA2A0_pxt { get; set; }
            //public double? I0_t_fromA2A1_pxt { get; set; }

        }
        public class SingleIsotopeCluster
        {
            float[] fIsotopeCluster = new float[6];
        }
        public struct TheoreticalI0Value
        {
            public string peptideseq;
            public double time;
            public double value;
            public double? charge;

            public TheoreticalI0Value(string peptideSeq, double t, double val, double charge)
            {
                peptideseq = peptideSeq;
                time = t;
                value = val;
                this.charge = charge;
            }

        }
    }
}

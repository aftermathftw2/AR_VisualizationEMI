using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Text;
using System.Threading;
using UnityEngine;

public class DataHandler : MonoBehaviour
{
    public GameObject dataPointMarker;
    public HandTracking handTracker;

    private PicoScopeHandler psHandler;

    [Range(0, 1)]
    public float measurementRadius = 0.07f;

    public List<MeasurementPoint> measurementPoints;
    public List<MeasurementPoint> visualizedMeasurementPoints;

    int frameLimiterIndex;
    int measurementIndex;
    int visualizationIndex;

    System.Tuple<double, double, TransformData> latestMeasurementData;
    private bool dataThreadRunning;
    private float radiusOfValidData = 0.035f;

    private string outputFile = "DataOutput.txt";

    // Awake is called upon initialization of the script
    void Awake()
    {
        frameLimiterIndex = 0;
        measurementIndex = -1;
        visualizationIndex = -1;

        measurementPoints = new List<MeasurementPoint>();
        visualizedMeasurementPoints = new List<MeasurementPoint>();

        //booting the picoscope
        Tuple<short, bool> setupStuff;
        try
        {
            setupStuff = PicoScopeHandler.PerformSetup();
            psHandler = new PicoScopeHandler(setupStuff.Item1, setupStuff.Item2);
            //start data collection thread
            StartCoroutine(waitForSecondsCoroutine(5));
            psHandler.StartThread();
            psHandler.usePingPong = false;
        }
        catch (Exception e)
        {
            Debug.Log(String.Format("unable to open PicoScope {0}", e.Message));
            Debug.Log("reverting to pingpong value");
            psHandler = new PicoScopeHandler(true);
        }
        StartCoroutine(waitForSecondsCoroutine(5));
        StartDataThread();
    }

    IEnumerator waitForSecondsCoroutine(float seconds)
    {
        //yield on a new YieldInstruction that waits for 5 seconds.
        yield return new WaitForSeconds(seconds);
    }

    private void OnDestroy()
    {
        dataThreadRunning = false;
        psHandler.StopPico(psHandler._handle);
    }

    void OnApplicationQuit()
    {
        dataThreadRunning = false;
        psHandler.StopPico(psHandler._handle);
        SaveDataFile();
    }

    void DataThreadRun()
    {
        while (dataThreadRunning)
        {
            while (!psHandler.doneCopying)
            {
                Thread.Sleep(10);
            }
            TransformData handTransFormData = new TransformData();
            if (handTracker.HandDetected)
            {
                //var handTransform = handTracker.GetTransform();
                var handTransform = handTracker.latestHandTransform;
                if (handTransform != null)
                {
                    handTransFormData = handTransform;
                }
            }
            System.Diagnostics.Stopwatch stopWatch = new System.Diagnostics.Stopwatch();
            stopWatch.Start();
            System.Tuple<double, double> measurementData = psHandler.GetData();
            latestMeasurementData = Tuple.Create(measurementData.Item1, measurementData.Item2, handTransFormData);
            stopWatch.Stop();
            //Debug.Log(String.Format("post processing done in {0} Milliseconds", stopWatch.ElapsedMilliseconds));
        }
    }

    public void StartDataThread()
    {
        Thread dataProcessThread = new Thread(DataThreadRun)
        {
            Name = "dataProcessThread"
        };
        dataThreadRunning = true;
        dataProcessThread.Start();
    }

    private void AddMeasurement()
    {
        
            measurementIndex++;
        TransformData handTransFormData = new TransformData(handTracker.GetTransform());

            /* This causes crashes because current200kpackage gets wiped every second. 
             * if the initialization of 200kpackage is done in the constructor it might be fine.
            while (!IsInRadius(latestMeasurementData.Item3.LocalPosition, handTransFormData.LocalPosition, radiusOfValidData))
            {
                StartCoroutine(waitForSecondsCoroutine(0.25f));
                Debug.Log("Waiting for new measurement point");
                //wait until new data is available which is closer to the actual hand
                //this does mean it can get stuck if you keep moving your hand....
                //and means massive lag spike...
            }
            */
            System.Tuple<double, double> measurementData = Tuple.Create(latestMeasurementData.Item1, latestMeasurementData.Item2);
            measurementPoints.Add(new MeasurementPoint() { transformData = handTransFormData, MeasurementData = measurementData.Item1, FrequencyOfMax = measurementData.Item2, MeasurementIndex = measurementIndex });
        
    }

    private void Visualize(int indexNewMeasurement)
    {
        //For performance this shouldnt loop over all measurements, only ones that are nearby currentMeasurementPoint
        //  how to split it up in chunks? needs to first check which chunk it is in before checking against all other points? what about chunk boundaries? check nearby chunks as well?
        // for now skip this.
        //should be split up in 2 (probably more) parts.
        //  One checking the just made measurementPoint if it is close enough to another point that it needs to overwrite.
        //  and one checking the measurementPoints and interpolating between them

        MeasurementPoint newMeasurementPoint = measurementPoints[indexNewMeasurement];
        bool inRadius = false;

        foreach (MeasurementPoint aPoint in visualizedMeasurementPoints)
        {
            //check first against visualizedMeasurementPoints if it overwrites previous visualizedPoints.
            if (IsInRadius(newMeasurementPoint.transformData.LocalPosition, aPoint.transformData.LocalPosition, measurementRadius))
            {
                inRadius = true;
                //it is in radius. if it is higher it should overwrite and update visuals
                if (newMeasurementPoint.MeasurementData > aPoint.MeasurementData)
                {
                    aPoint.MeasurementData = newMeasurementPoint.MeasurementData;
                    aPoint.FrequencyOfMax = newMeasurementPoint.FrequencyOfMax;
                    
                    aPoint.Visualization.GetComponent<Renderer>().material.color = determineLerpColour(aPoint);
                    //Debug.Log(string.Format("{0}, {1}, higher value", aPoint, inRadius));
                }
                //it is in radius. but it is lower. keep original visualizationPoint
                else
                {
                    //Debug.Log(string.Format("{0}, {1}, lower value", aPoint, inRadius));
                }
            }
        }
        if (!inRadius)
        {
            Debug.Log(string.Format("{0}, {1}, new point", newMeasurementPoint, inRadius));
            visualizationIndex++;
            visualizedMeasurementPoints.Add(new MeasurementPoint() {transformData = newMeasurementPoint.transformData, MeasurementData = newMeasurementPoint.MeasurementData, MeasurementIndex = visualizationIndex, Visualization = Instantiate(dataPointMarker, newMeasurementPoint.transformData.LocalPosition, Quaternion.identity)});

            visualizedMeasurementPoints[visualizationIndex].Visualization.GetComponent<Renderer>().material.color = determineLerpColour(newMeasurementPoint);
        }
        //now either the point is discarded, older points are updated or a new point is made.
        //now interpolation between existing points can happen.
        
    }

    private Color determineLerpColour(MeasurementPoint aPoint)
    {
        double standardAtFreq = 1000;

        for (int i = 0; i < psHandler.post.standard.Item1.Count; i++)
        {
            if (aPoint.FrequencyOfMax == (double)psHandler.post.standard.Item1[i])
            {
                standardAtFreq = psHandler.post.standard.Item2[i];
                i = psHandler.post.standard.Item1.Count - 1;
            }
        }
        //Debug.Log(String.Format("standard: {0} at freq: {1}", standardAtFreq, aPoint.FrequencyOfMax));
        Color lerpColour = Color.Lerp(Color.red, Color.yellow, (float)((aPoint.MeasurementData) / standardAtFreq));
        if (standardAtFreq - aPoint.MeasurementData > 20)
        {
            //lerp between black and red if the distance from the limitline is bigger than 20
            lerpColour = Color.Lerp(Color.black, Color.red, (float)((aPoint.MeasurementData - standardAtFreq + 50) / 30));
            lerpColour.a = 0.3f;
        }
        else if (standardAtFreq - aPoint.MeasurementData > 0 && standardAtFreq - aPoint.MeasurementData < 20)
        {
            //lerp between red and yellow if the distance from the limitline is bigger than 0 but smaller than 20
            lerpColour = Color.Lerp(Color.red, Color.yellow, (float)((aPoint.MeasurementData - standardAtFreq + 20) / 20));
            lerpColour.a = 0.6f;
        }
        else if (standardAtFreq - aPoint.MeasurementData < 0)
        {
            //lerp between yellow and white if it is above the limitline
            lerpColour = Color.Lerp(Color.yellow, Color.white, (float)((aPoint.MeasurementData - standardAtFreq + 0) / 10));
            lerpColour.a = 0.9f;
        }
        
        return lerpColour;
    }

    /// <summary>
    /// Checks if point 2 is in radius of point 1. 
    /// </summary>
    /// <param name="point1"></param>
    /// <param name="point2"></param>
    /// <param name="radius"></param>
    /// <returns>True if it is in range, false if not </returns>
    private bool IsInRadius(Vector3 point1, Vector3 point2, float radius)
    {
        float R = Mathf.Sqrt(Mathf.Pow(point1.x - point2.x, 2) + Mathf.Pow(point1.y - point2.y, 2) + Mathf.Pow(point1.z - point2.z, 2));
        if (R < radius)
        {
            return true;
        }
        return false;
    }

    private void SaveDataFile()
    {
        //loop over all visualized measurement points
        //save their x, y, z, value and frequency
        var sb = new StringBuilder();
        
        foreach (MeasurementPoint mPoint in visualizedMeasurementPoints)
        {
            sb.AppendLine(String.Format("{0}, {1}, {2}, {3}, {4}", 
                mPoint.transformData.LocalPosition.x,
                mPoint.transformData.LocalPosition.y,
                mPoint.transformData.LocalPosition.z,
                mPoint.MeasurementData,
                mPoint.FrequencyOfMax));
        }

        using (TextWriter writer = new StreamWriter(outputFile, false))
        {
            writer.Write(sb.ToString());
            writer.Close();
        }
    }

    //update gets called every frame
    private void Update()
    {
        if (handTracker.HandDetected)
        {
            frameLimiterIndex++;
            if (frameLimiterIndex % 10 == 0)
            {
                if (latestMeasurementData != null)
                {
                    AddMeasurement();
                    Visualize(measurementIndex);
                }
            }

        }
    }
}






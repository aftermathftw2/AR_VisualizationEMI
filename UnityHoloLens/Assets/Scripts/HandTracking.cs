using System.Collections;
using System.Collections.Generic;
using UnityEngine;

using System.Linq;
using UnityEngine.XR.WSA.Input;
using System;

public class HandTracking : MonoBehaviour
{

    public GameObject sphereTrackingMarker;

    private HashSet<uint> trackedHands = new HashSet<uint>();
    private Dictionary<uint, GameObject> trackingObject = new Dictionary<uint, GameObject>();
    private GestureRecognizer gestureRecognizer;
    private uint activeId;
    private bool holding;

    public TransformData latestHandTransform;

    public bool HandDetected
    {
        get { return trackedHands.Count > 0; }
    }

    // Awake is called on initializing of the script
    void Awake()
    {
        InteractionManager.InteractionSourceDetected += InteractionManager_InteractionSourceDetected;
        InteractionManager.InteractionSourceUpdated += InteractionManager_InteractionSourceUpdated;
        InteractionManager.InteractionSourceLost += InteractionManager_InteractionSourceLost;

        gestureRecognizer = new GestureRecognizer();
#pragma warning disable CS0618 // Type or member is obsolete
        gestureRecognizer.SetRecognizableGestures(GestureSettings.Tap | GestureSettings.Hold);

        gestureRecognizer.Tapped += GestureRecognizerTapped;
        gestureRecognizer.HoldStarted += GestureRecognizer_HoldStarted;
        gestureRecognizer.HoldCompleted += GestureRecognizer_HoldCompleted;
        gestureRecognizer.HoldCanceled += GestureRecognizer_HoldCanceled;
        gestureRecognizer.StartCapturingGestures();
        latestHandTransform = new TransformData();
    }

    private void InteractionManager_InteractionSourceDetected(InteractionSourceDetectedEventArgs args)
    {
        uint id = args.state.source.id;
        // Check to see that the source is a hand.
        if (args.state.source.kind != InteractionSourceKind.Hand)
        {
            return;
        }
        trackedHands.Add(id);
        activeId = id;

        var obj = Instantiate(sphereTrackingMarker);
        Vector3 pos;

        if (args.state.sourcePose.TryGetPosition(out pos))
        {
            obj.transform.position = pos;
            latestHandTransform.LocalPosition = pos;
        }

        trackingObject.Add(id, obj);
    }

    private void InteractionManager_InteractionSourceUpdated(InteractionSourceUpdatedEventArgs args)
    {
        uint id = args.state.source.id;
        Vector3 pos;
        Quaternion rot;

        if (args.state.source.kind == InteractionSourceKind.Hand)
        {
            if (trackingObject.ContainsKey(id))
            {
                if (args.state.sourcePose.TryGetPosition(out pos))
                {
                    trackingObject[id].transform.position = pos;
                    latestHandTransform.LocalPosition = pos;
                }

                if (args.state.sourcePose.TryGetRotation(out rot))
                {
                    trackingObject[id].transform.rotation = rot;
                    latestHandTransform.LocalEulerRotation = rot.eulerAngles;
                }
            }
        }
    }

    private void InteractionManager_InteractionSourceLost(InteractionSourceLostEventArgs args)
    {
        uint id = args.state.source.id;
        // Check to see that the source is a hand.
        if (args.state.source.kind != InteractionSourceKind.Hand)
        {
            return;
        }

        if (trackedHands.Contains(id))
        {
            trackedHands.Remove(id);
        }

        if (trackingObject.ContainsKey(id))
        {
            var obj = trackingObject[id];
            trackingObject.Remove(id);
            Destroy(obj);
        }
        if (trackedHands.Count > 0)
        {
            activeId = trackedHands.First();
        }
    }

    private void GestureRecognizerTapped(TappedEventArgs args)
    {
        uint id = args.source.id;
        if (trackingObject.ContainsKey(activeId))
        {
            Debug.Log("tapped");
        }
    }

    private void GestureRecognizer_HoldStarted(HoldStartedEventArgs args)
    {
        uint id = args.source.id;
        if (trackingObject.ContainsKey(activeId))
        {
            holding = true;
        }
    }

    private void GestureRecognizer_HoldCompleted(HoldCompletedEventArgs args)
    {
        uint id = args.source.id;
        if (trackingObject.ContainsKey(activeId))
        {
            holding = false;
        }
    }

    private void GestureRecognizer_HoldCanceled(HoldCanceledEventArgs args)
    {
        uint id = args.source.id;
        if (trackingObject.ContainsKey(activeId))
        {
            holding = false;
        }
    }

    public bool IsHolding()
    {
        return holding;
    }
    
    public Transform GetTransform()
    {
        if (trackingObject.ContainsKey(activeId))
        {
            return trackingObject[activeId].transform;
        } else
        {
            Debug.Log("No tracked Hand to return transform of");
            Debug.Log(trackingObject.Count);
            return null;
        }
    }


    private void OnDestroy()
    {
        InteractionManager.InteractionSourceDetected -= InteractionManager_InteractionSourceDetected;
        InteractionManager.InteractionSourceUpdated -= InteractionManager_InteractionSourceUpdated;
        InteractionManager.InteractionSourceLost -= InteractionManager_InteractionSourceLost;

        gestureRecognizer.Tapped -= GestureRecognizerTapped;
        gestureRecognizer.HoldStarted -= GestureRecognizer_HoldStarted;
        gestureRecognizer.HoldCompleted -= GestureRecognizer_HoldCompleted;
        gestureRecognizer.HoldCanceled -= GestureRecognizer_HoldCanceled;
        gestureRecognizer.StopCapturingGestures();
    }
}

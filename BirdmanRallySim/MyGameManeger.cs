using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class MyGameManeger : MonoBehaviour
{
    public static MyGameManeger instance = null;
    
    [System.NonSerialized] public bool Landing = false;
    [System.NonSerialized] public bool HUDActive = true;
    [System.NonSerialized] public bool HorizontalLineActive = false;
    [System.NonSerialized] public bool SettingActive = false;
    [System.NonSerialized] public bool CameraSwitch = true; // true:FPS false:TPS
    [System.NonSerialized] public bool SettingChanged = false;
    [System.NonSerialized] public bool MousePitchControl = false;
    [System.NonSerialized] public float MouseSensitivity = 1.000f; // Magnitude of Gust [m/s]
    [System.NonSerialized] public float GustMag = 0.000f; // Magnitude of Gust [m/s]
    [System.NonSerialized] public float GustDirection = 0.000f; // Direction of Gust [deg]: -180~180
    [System.NonSerialized] public float Airspeed_TO = 5.000f; // Airspeed at take-off [m/s]
    [System.NonSerialized] public float alpha_TO = 0.000f; // Angle of attack at take-off [deg]
    [System.NonSerialized] public string PlaneName = "QX-20";
    [System.NonSerialized] public string FlightMode = "BirdmanRally";
    [System.NonSerialized] public GameObject Plane = null;
    [System.NonSerialized] public Vector3 PlatformPosition = new Vector3(0f,10.5f,0f);

    // Start is called before the first frame update
    void Awake()
    {
        if(instance == null)
        {
            instance = this;
            DontDestroyOnLoad(this.gameObject); 
        }
        else
        {
            Destroy(this.gameObject);
        }
    }
}

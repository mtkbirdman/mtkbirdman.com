using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class MyGameManeger : MonoBehaviour
{
    public static MyGameManeger instance = null;
    public bool Landing = false;
    public bool HUDActive = true;
    public bool HorizontalLineActive = false;
    public bool SettingActive = false;
    public bool CameraSwitch = false;
    public bool SettingChanged = false;
    public float GustMag = 0.000f; // Magnitude of Gust [m/s]
    public float GustDirection = 0.000f; // Direction of Gust [deg]: -180~180
    public float Airspeed_TO = 0.000f; // Airspeed at take-off [m/s]
    public float alpha_TO = 0.000f; // Angle of attack at take-off [deg]
    public string PlaneName = "QX-20";
    public GameObject Plane = null;

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

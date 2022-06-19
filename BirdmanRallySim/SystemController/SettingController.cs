using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;
using UnityEngine.SceneManagement;
using System;

public class SettingController : MonoBehaviour
{
    private GameObject Setting;
    private Slider TakeoffVelocitySlider;

    // Start is called before the first frame update
    void Start()
    {
        Setting = GameObject.Find("Setting");
        TakeoffVelocitySlider = GameObject.Find("TakeoffVelocitySlider").GetComponent<Slider>();

        MyGameManeger.instance.Airspeed_TO = TakeoffVelocitySlider.value*0.1f;
        
        MyGameManeger.instance.SettingActive = false;
        Setting.SetActive(MyGameManeger.instance.SettingActive);
    }

    // Update is called once per frame
    void Update()
    {
        if(Input.GetKeyDown("tab")){
            MyGameManeger.instance.SettingActive = !MyGameManeger.instance.SettingActive;
            Setting.SetActive(MyGameManeger.instance.SettingActive);
            Time.timeScale=(float)Convert.ToInt32(!MyGameManeger.instance.SettingActive & !MyGameManeger.instance.Landing);
        }
        if(Input.GetKeyDown("c")){
            MyGameManeger.instance.MousePitchControl = !MyGameManeger.instance.MousePitchControl;
        }
    }
}

using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class SettingController : MonoBehaviour
{
    private GameObject SETTING; // Canvas

    // Start is called before the first frame update
    void Start()
    {
        SETTING = GameObject.Find("SETTING");

        SETTING.SetActive(false);
        MyGameManeger.instance.SettingActive = false;
    }

    // Update is called once per frame
    void Update()
    {
        if(Input.GetKeyDown("tab")){
            if(MyGameManeger.instance.SettingActive){ // Close
                SETTING.SetActive(false);
                MyGameManeger.instance.SettingActive = false;
                if(!MyGameManeger.instance.Landing){Time.timeScale=1f;}
            }
            else{ // Open
                SETTING.SetActive(true);
                MyGameManeger.instance.SettingActive = true;
                Time.timeScale=0f;
            }
        }
    }
}

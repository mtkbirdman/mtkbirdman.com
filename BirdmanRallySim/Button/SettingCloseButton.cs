using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;
using UnityEngine.SceneManagement;
using System;

public class SettingCloseButton : MonoBehaviour
{
    private GameObject Setting;
    private bool firstPush = false;

    // Start is called before the first frame update
    void Start()
    {
        Setting = GameObject.Find("Setting");
    }

    public void OnClick()
    {
        if (!firstPush)
        {
            MyGameManeger.instance.SettingActive = !MyGameManeger.instance.SettingActive;
            Setting.SetActive(MyGameManeger.instance.SettingActive);
            Time.timeScale=(float)Convert.ToInt32(!MyGameManeger.instance.SettingActive & !MyGameManeger.instance.Landing);
        }
    }
}

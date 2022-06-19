using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.SceneManagement;
using System;

public class Landing : MonoBehaviour
{
    private GameObject Result;
    
    // Start is called before the first frame update
    void Start()
    {
        Result = GameObject.Find("Result");

        Result.SetActive(false);
        MyGameManeger.instance.Landing = false;
    }

    // Update is called once per frame
    void OnCollisionEnter(Collision collision)
    {
        MyGameManeger.instance.Landing = true;
        Time.timeScale=(float)Convert.ToInt32(!MyGameManeger.instance.SettingActive & !MyGameManeger.instance.Landing);
    }

    // Update is called once per frame
    void Update()
    {
        Result.SetActive(!MyGameManeger.instance.SettingActive & MyGameManeger.instance.Landing);
    }
}

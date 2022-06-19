using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System;

public class HUDController : MonoBehaviour
{
    //private GameObject HUD;
    private GameObject HorizontalLine;
    private GameObject HUDCanvas;

    // Start is called before the first frame update
    void Start()
    {
        HUDCanvas = GameObject.Find("HUD");
        HorizontalLine = GameObject.Find("HUD").transform.Find("HorizontalLine").gameObject;
        
        HUDCanvas.SetActive(MyGameManeger.instance.HUDActive);
        HorizontalLine.SetActive(MyGameManeger.instance.HorizontalLineActive);
    }

    // Update is called once per frame
    void Update()
    {
        if(Input.GetKeyDown("h")){
            MyGameManeger.instance.HUDActive = !MyGameManeger.instance.HUDActive;
            HUDCanvas.SetActive(MyGameManeger.instance.HUDActive);
        }

        if(MyGameManeger.instance.CameraSwitch && Input.GetKeyDown("l")){
            HorizontalLine.SetActive(!MyGameManeger.instance.HorizontalLineActive);
            MyGameManeger.instance.HorizontalLineActive = !MyGameManeger.instance.HorizontalLineActive;
        }
    }
}

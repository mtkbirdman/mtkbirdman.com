using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class MyCameraController : MonoBehaviour
{
    private GameObject HorizontalLine;
    private Camera FPSCamera; // Camera
    private Camera TPSCamera; // Camera
    private Canvas HUDCanvas;
    private Canvas SETTINGCanvas;

    // Start is called before the first frame update
    void Start()
    {
        FPSCamera = MyGameManeger.instance.Plane.transform.Find("FPSCamera").gameObject.GetComponent<Camera>();
        TPSCamera = MyGameManeger.instance.Plane.transform.Find("TPSCamera").gameObject.GetComponent<Camera>();
        HUDCanvas = GameObject.Find("HUD").GetComponent<Canvas>();
        SETTINGCanvas = GameObject.Find("SETTING").GetComponent<Canvas>();
        HorizontalLine = GameObject.Find("HUD").transform.Find("HorizontalLine").gameObject;

        if(MyGameManeger.instance.CameraSwitch){
            MyGameManeger.instance.CameraSwitch = false;
            SwitchCamera(); // TPS to FPS
        }else{
            MyGameManeger.instance.CameraSwitch = true;
            SwitchCamera(); // TPS to TPS
        }
    }

    // Update is called once per frame
    void Update()
    {
        if(Input.GetKeyDown("space")){SwitchCamera();}
    }

    void SwitchCamera()
    {
        if(MyGameManeger.instance.CameraSwitch){ 
            // FPS to TPS
            TPSCamera.enabled = true;
            FPSCamera.enabled = false;
            HUDCanvas.worldCamera = TPSCamera;
            SETTINGCanvas.worldCamera = TPSCamera;
            MyGameManeger.instance.CameraSwitch = false;
            HorizontalLine.SetActive(false);
            MyGameManeger.instance.HorizontalLineActive = false;
        }
        else{ 
            // TPS to FPS
            TPSCamera.enabled = false;
            FPSCamera.enabled = true;
            HUDCanvas.worldCamera = FPSCamera;
            SETTINGCanvas.worldCamera = FPSCamera;
            MyGameManeger.instance.CameraSwitch = true;
            if(MyGameManeger.instance.HorizontalLineActive){
                HorizontalLine.SetActive(true);
                MyGameManeger.instance.HorizontalLineActive = true;
            }
        }
    }
}

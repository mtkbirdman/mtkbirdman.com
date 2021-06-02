using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class HUDController : MonoBehaviour
{
    //private GameObject HUD;
    private GameObject HorizontalLine;
    private Canvas HUDCanvas;

    // Start is called before the first frame update
    void Start()
    {
        //HUD = GameObject.Find("HUD");
        HUDCanvas = GameObject.Find("HUD").GetComponent<Canvas>();
        HorizontalLine = GameObject.Find("HUD").transform.Find("HorizontalLine").gameObject;

        if(MyGameManeger.instance.HUDActive){
            HUDCanvas.planeDistance = 0.11f;
            MyGameManeger.instance.HUDActive = true;
            if(MyGameManeger.instance.HorizontalLineActive){
                HorizontalLine.SetActive(true);
                MyGameManeger.instance.HorizontalLineActive = true;
            }else{
                HorizontalLine.SetActive(false);
                MyGameManeger.instance.HorizontalLineActive = false;
            }
        }else{
            HUDCanvas.planeDistance = 0f;
            MyGameManeger.instance.HUDActive = false;
        }
    }

    // Update is called once per frame
    void Update()
    {
        if(MyGameManeger.instance.HUDActive){
            if(Input.GetKeyDown("i")){
                HUDCanvas.planeDistance = 0f;
                MyGameManeger.instance.HUDActive = false;
            }
            if(MyGameManeger.instance.CameraSwitch && Input.GetKeyDown("h")){
                if(MyGameManeger.instance.HorizontalLineActive){
                    HorizontalLine.SetActive(false);
                    MyGameManeger.instance.HorizontalLineActive = false;
                }else{
                    HorizontalLine.SetActive(true);
                    MyGameManeger.instance.HorizontalLineActive = true;
                }
            }
        }else{
            if(Input.GetKeyDown("i")){
                HUDCanvas.planeDistance = 0.11f;
                MyGameManeger.instance.HUDActive = true;
            }
        }
    }
}

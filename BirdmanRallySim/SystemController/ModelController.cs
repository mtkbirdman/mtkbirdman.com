using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class ModelController : MonoBehaviour
{
    private GameObject PlaneParent;
    
    // Start is called before the first frame update
    void Start()
    {
        PlaneParent = GameObject.Find("Plane");
        foreach(Transform item in PlaneParent.transform){
            if(item.gameObject.name != MyGameManeger.instance.PlaneName){
                item.gameObject.SetActive(false);
            }
            if(item.gameObject.name == "WaterProDaytime"){
                item.gameObject.SetActive(true);
            }
        }
        MyGameManeger.instance.Plane = GameObject.Find(MyGameManeger.instance.PlaneName);
        //Debug.Log(MyGameManeger.instance.Plane);
    }
}

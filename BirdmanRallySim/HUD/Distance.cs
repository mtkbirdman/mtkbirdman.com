using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class Distance : MonoBehaviour
{
    private Text scoreText;
    private Rigidbody PlaneRigidbody;
    
    
    // Start is called before the first frame update
    void Start()
    {
        scoreText = this.GetComponent<Text>();
        PlaneRigidbody = MyGameManeger.instance.Plane.GetComponent<Rigidbody>();
    }

    // Update is called once per frame
    void Update()
    {
        float Distance = (PlaneRigidbody.position-MyGameManeger.instance.PlatformPosition).magnitude;
        if(MyGameManeger.instance.FlightMode=="BirdmanRally"){Distance-=10f;}

        scoreText.text = "\n" + Distance.ToString("0.000");
    }
}

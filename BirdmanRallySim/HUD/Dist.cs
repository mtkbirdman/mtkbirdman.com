using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class Dist : MonoBehaviour
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
        scoreText.text = "\n" + PlaneRigidbody.position.magnitude.ToString("0.000");
    }
}

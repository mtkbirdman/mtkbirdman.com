using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class Rotation : MonoBehaviour
{
    private Text scoreText;

    // Start is called before the first frame update
    void Start()
    {
        scoreText = this.GetComponent<Text>();
    }

    void Update()
    {
        // Calculate rotation
        float q1 = MyGameManeger.instance.Plane.transform.rotation.x;
        float q2 = -MyGameManeger.instance.Plane.transform.rotation.y;
        float q3 = -MyGameManeger.instance.Plane.transform.rotation.z;
        float q4 = MyGameManeger.instance.Plane.transform.rotation.w;
        float C11 = q1*q1-q2*q2-q3*q3+q4*q4;
        float C22 = -q1*q1+q2*q2-q3*q3+q4*q4;
        float C12 = 2f*(q1*q2+q3*q4);
        float C13 = 2f*(q1*q3-q2*q4);
        float C32 = 2f*(q2*q3-q1*q4);
        float phi = -Mathf.Atan(-C32/C22)*Mathf.Rad2Deg;
        float theta = -Mathf.Asin(C12)*Mathf.Rad2Deg; 
        float psi = -Mathf.Atan(-C13/C11)*Mathf.Rad2Deg;

        scoreText.text = 
            "\n\n" + phi.ToString("0.000") 
            + "\r\n" + theta.ToString("0.000") 
            + "\r\n" + psi.ToString("0.000");
    }
    
}

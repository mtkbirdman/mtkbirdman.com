using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class AerodynamicCalculator : MonoBehaviour
{
    // public
    [System.NonSerialized] public float Airspeed = 0.000f; // Airspeed [m/s]
    [System.NonSerialized] public float alpha = 0.000f; // Angle of attack [deg]
    [System.NonSerialized] public float beta = 0.000f; // Side slip angle [deg]
    [System.NonSerialized] public float de = 0.000f; // Elevator angle [deg]
    [System.NonSerialized] public float dr = 0.000f; // Rudder angle [deg]
    [System.NonSerialized] public float dh = 0.000f; // Movement of c.g. [-]
    [System.NonSerialized] public float LocalGustMag = 0.000f; // Magnitude of local gust [m/s]
    [System.NonSerialized] public float LocalGustDirection = 0.000f; // Magnitude of local gust [m/s]
    // Phisics
    static private float rho = 1.164f;
    static private float hE0 = 10.500f; // Altitude at Take-off [m]
    // At Cruise without Ground Effect
    static private float Airspeed0; // Magnitude of ground speed [m/s]
    static private float alpha0; // Angle of attack [deg]
    static private float CDp0; // Parasitic drag [-]
    static private float Cmw0; // Pitching momentum [-]
    static private float CLMAX; // Lift Coefficient [-]
    static private float CL0 = 0.000f; // Lift Coefficient [-]
    static private float CLw0 = 0.000f; // Lift Coefficient [-]
    static private float CLt0 = 0.000f; // Tail Coefficient [-]
    static private float epsilon0 = 0.000f; // Downwash
    // Plane
    static bool Downwash; // Conventional Tail: True, T-Tail: False
    static private float CL = 0.000f; // Lift Coefficient [-]
    static private float CD = 0.000f; // Drag Coefficient [-]
    static private float Cx = 0.000f; // X Force Coefficient [-]
    static private float Cy = 0.000f; // Y Force Coefficient [-]
    static private float Cz = 0.000f; // Z Force Coefficient [-]
    static private float Cl = 0.000f; // Rolling momentum [-]
    static private float Cm = 0.000f; // Pitching momentum [-]
    static private float Cn = 0.000f; // Yawing momentum [-]
    static private float dh0 = 0.000f; // Initial Mouse Position
    // Wing
    static private float Sw; // Wing area of wing [m^2]
    static private float bw; // Wing span [m]
    static private float cMAC; // Mean aerodynamic chord [m]
    static public float aw; // Wing Lift Slope [1/deg]
    static private float hw; // Length between Wing a.c. and c.g. [-]
    static private float AR; // Aspect Ratio [-]
    static private float ew; // Wing efficiency [-]
    static private float CLw = 0.000f; // Lift Coefficient [-]
    // Tail
    static private float St; // Wing area of tail [m^2]
    static private float at; // Tail Lift Slope [1/deg]
    static private float lt; // Length between Tail a.c. and c.g. [m]
    static private float VH; // Tail Volume [-]
    static private float deMAX; // Maximum elevator angle [deg]
    static private float tau; // Control surface angle of attack effectiveness [-]
    static private float CLt = 0.000f; // Lift Coefficient [-]
    // Fin
    static private float drMAX; // Maximum rudder angle
    // Ground Effect
    static private float CGEMIN; // Minimum Ground Effect Coefficient [-]
    static private float CGE = 0f; // Ground Effect Coefficient: CDiGE/CDi [-]
    // Stability derivatives
    static private float Cyb; // [1/deg]
    static private float Cyp; // [1/rad]
    static private float Cyr; // [1/rad]
    static private float Cydr; // [1/deg]
    static private float Cnb; // [1/deg]
    static private float Cnp; // [1/rad]
    static private float Cnr; // [1/rad]
    static private float Cndr; // [1/deg]
    static private float Clb; // [1/deg]
    static private float Clp; // [1/rad]
    static private float Clr; // [1/rad]
    static private float Cldr; // [1/deg]
    // Gust
    static private Vector3 Gust = Vector3.zero; // Gust [m/s]
    // Rotation
    static private float phi; // [deg]
    static private float theta;  // [deg]
    static private float psi; // [deg]

    private Rigidbody PlaneRigidbody;

    // Start is called before the first frame update
    void Start()
    {        
        // Get rigidbody component
        PlaneRigidbody = this.GetComponent<Rigidbody>();

        // Input Specifications
        InputSpecifications();
        
        // Set take-off speed
        if(MyGameManeger.instance.FlightMode=="BirdmanRally"){
            //MyGameManeger.instance.Airspeed_TO = 5.0f; // Airspeed at take-off [m/s]
            PlaneRigidbody.velocity = Vector3.zero;
        }else if(MyGameManeger.instance.FlightMode=="TestFlight"){ // 
            PlaneRigidbody.velocity = new Vector3(
                Airspeed0*Mathf.Cos(Mathf.Deg2Rad*alpha0),
                -Airspeed0*Mathf.Sin(Mathf.Deg2Rad*alpha0),
                0f
            );
        }

        // Calculate CL at cluise
        CL0 = (PlaneRigidbody.mass*Physics.gravity.magnitude)/(0.5f*rho*Airspeed0*Airspeed0*Sw);
        CLt0 = (Cmw0+CL0*hw)/(VH+(St/Sw)*hw);
        CLw0 = CL0-(St/Sw)*CLt0;
        if(Downwash){epsilon0 = (CL0/(Mathf.PI*ew*AR))*Mathf.Rad2Deg;}

        dh0 = Screen.height/2f; // Initial Mouse Position
    }
    
    void FixedUpdate()
    {
        
        // Velocity and AngularVelocity
        float u = transform.InverseTransformDirection(PlaneRigidbody.velocity).x;
        float v = -transform.InverseTransformDirection(PlaneRigidbody.velocity).z;
        float w = -transform.InverseTransformDirection(PlaneRigidbody.velocity).y;
        float p = -transform.InverseTransformDirection(PlaneRigidbody.angularVelocity).x*Mathf.Rad2Deg;
        float q = transform.InverseTransformDirection(PlaneRigidbody.angularVelocity).z*Mathf.Rad2Deg;
        float r = transform.InverseTransformDirection(PlaneRigidbody.angularVelocity).y*Mathf.Rad2Deg;
        float hE = PlaneRigidbody.position.y;

        // Force and Momentum
        Vector3 AerodynamicForce = Vector3.zero;
        Vector3 AerodynamicMomentum = Vector3.zero;
        Vector3 TakeoffForce = Vector3.zero;

        // Hoerner and Borst (Modified)
        CGE = (CGEMIN+33f*Mathf.Pow((hE/bw),1.5f))/(1f+33f*Mathf.Pow((hE/bw),1.5f));

        // Get control surface angles
        de = 0.000f;
        dr = 0.000f;
        if (MyGameManeger.instance.MousePitchControl){
            dh = -(Input.mousePosition.y-dh0)*0.0001f*MyGameManeger.instance.MouseSensitivity;
        }
        //Debug.Log(dh);

        de = Input.GetAxisRaw("Vertical")*deMAX;
        dr = -Input.GetAxisRaw("Horizontal")*drMAX;

        if(Input.GetMouseButton(0)){dr = drMAX;}
        else if(Input.GetMouseButton(1)){dr = -drMAX;}

        //Debug.Log(dh);

        // Gust
        LocalGustMag = MyGameManeger.instance.GustMag*Mathf.Pow((hE/hE0),1f/7f);
        Gust = Quaternion.AngleAxis(MyGameManeger.instance.GustDirection,Vector3.up)*(Vector3.right*LocalGustMag);
        Vector3 LocalGust = this.transform.InverseTransformDirection(Gust);
        float ug = LocalGust.x + 1e-10f;
        float vg = -LocalGust.z;
        float wg = -LocalGust.y;
        if(ug>0){LocalGustDirection = Mathf.Atan(vg/(ug+1e-10f))*Mathf.Rad2Deg;}
        else{LocalGustDirection = Mathf.Atan(vg/(ug+1e-10f))*Mathf.Rad2Deg+vg/Mathf.Abs((vg+1e-10f))*180;}

        // Calculate angles
        Airspeed = Mathf.Sqrt((u+ug)*(u+ug)+(v+vg)*(v+vg)+(w+wg)*(w+wg));
        alpha = Mathf.Atan((w+wg)/(u+ug))*Mathf.Rad2Deg;
        beta = Mathf.Atan((v+vg)/Airspeed)*Mathf.Rad2Deg;

        // Wing and Tail
        CLw = CLw0+aw*(alpha-alpha0);
        CLt = CLt0+at*((alpha-alpha0)+(1f-CGE*(CLw/CLw0))*epsilon0+de*tau+((lt-dh*cMAC)/Airspeed)*q);
        if(Mathf.Abs(CLw)>CLMAX){CLw = (CLw/Mathf.Abs(CLw))*CLMAX;} // Stall
        if(Mathf.Abs(CLt)>CLMAX){CLt = (CLt/Mathf.Abs(CLt))*CLMAX;} // Stall

        // Lift and Drag
        CL = CLw+(St/Sw)*CLt; // CL        
        CD = CDp0*(1f+Mathf.Abs(Mathf.Pow((alpha/9f),3f)))+((CL*CL)/(Mathf.PI*ew*AR))*CGE; // CD

        // Force
        Cx = CL*Mathf.Sin(Mathf.Deg2Rad*alpha)-CD*Mathf.Cos(Mathf.Deg2Rad*alpha); // Cx       
        Cy = Cyb*beta+Cyp*(1f/Mathf.Rad2Deg)*((p*bw)/(2f*Airspeed))+Cyr*(1f/Mathf.Rad2Deg)*((r*bw)/(2f*Airspeed))+Cydr*dr; // Cy       
        Cz = -CL*Mathf.Cos(Mathf.Deg2Rad*alpha)-CD*Mathf.Sin(Mathf.Deg2Rad*alpha); // Cz

        // Torque
        Cl = Clb*beta+Clp*(1f/Mathf.Rad2Deg)*((p*bw)/(2f*Airspeed))+Clr*(1f/Mathf.Rad2Deg)*((r*bw)/(2f*Airspeed))+Cldr*dr; // Cl        
        Cm = Cmw0+CLw*hw-VH*CLt+CL*dh; // Cm       
        Cn = Cnb*beta+Cnp*(1f/Mathf.Rad2Deg)*((p*bw)/(2f*Airspeed))+Cnr*(1f/Mathf.Rad2Deg)*((r*bw)/(2f*Airspeed))+Cndr*dr; // Cn

        AerodynamicForce.x = 0.5f*rho*Airspeed*Airspeed*Sw*Cx;
        AerodynamicForce.y = 0.5f*rho*Airspeed*Airspeed*Sw*(-Cz);
        AerodynamicForce.z = 0.5f*rho*Airspeed*Airspeed*Sw*(-Cy);

        AerodynamicMomentum.x = 0.5f*rho*Airspeed*Airspeed*Sw*bw*(-Cl);
        AerodynamicMomentum.y = 0.5f*rho*Airspeed*Airspeed*Sw*bw*Cn;
        AerodynamicMomentum.z = 0.5f*rho*Airspeed*Airspeed*Sw*cMAC*Cm;


        float Distance = (PlaneRigidbody.position-MyGameManeger.instance.PlatformPosition).magnitude-10f;
        if(MyGameManeger.instance.FlightMode=="BirdmanRally" && Distance<-0.5f){
            
            CalculateRotation();
            
            float W = PlaneRigidbody.mass*Physics.gravity.magnitude;
            float L = 0.5f*rho*Airspeed*Airspeed*Sw*(Cx*Mathf.Sin(Mathf.Deg2Rad*theta)-Cz*Mathf.Cos(Mathf.Deg2Rad*theta));
            float N = (W-L)*Mathf.Cos(Mathf.Deg2Rad*3.5f); // N=(W-L)*cos(3.5deg)
            float P = (PlaneRigidbody.mass*MyGameManeger.instance.Airspeed_TO*MyGameManeger.instance.Airspeed_TO)/(2f*10f); // P=m*Vto*Vto/2*L
            
            TakeoffForce.x = P;
            TakeoffForce.y = N*Mathf.Cos(Mathf.Deg2Rad*3.5f);
            TakeoffForce.z = 0f;
            
            AerodynamicForce.z = 0f;
            AerodynamicMomentum.x = 0f;
            AerodynamicMomentum.y = 0f;
        }
        //Debug.Log(Distance);

        PlaneRigidbody.AddRelativeForce(AerodynamicForce, ForceMode.Force);
        PlaneRigidbody.AddRelativeTorque(AerodynamicMomentum, ForceMode.Force);
        PlaneRigidbody.AddForce(TakeoffForce, ForceMode.Force);
    }

    void CalculateRotation()
    {
        float q1 = MyGameManeger.instance.Plane.transform.rotation.x;
        float q2 = -MyGameManeger.instance.Plane.transform.rotation.y;
        float q3 = -MyGameManeger.instance.Plane.transform.rotation.z;
        float q4 = MyGameManeger.instance.Plane.transform.rotation.w;
        float C11 = q1*q1-q2*q2-q3*q3+q4*q4;
        float C22 = -q1*q1+q2*q2-q3*q3+q4*q4;
        float C12 = 2f*(q1*q2+q3*q4);
        float C13 = 2f*(q1*q3-q2*q4);
        float C32 = 2f*(q2*q3-q1*q4);

        phi = -Mathf.Atan(-C32/C22)*Mathf.Rad2Deg;
        theta = -Mathf.Asin(C12)*Mathf.Rad2Deg; 
        psi = -Mathf.Atan(-C13/C11)*Mathf.Rad2Deg;
    }

    void InputSpecifications()
    {
        if(MyGameManeger.instance.PlaneName == "QX-18"){
            // Plane
            PlaneRigidbody.mass = 93.875f; // [kg]
            PlaneRigidbody.centerOfMass = new Vector3(0f,0.221f,0f); // [m]
            PlaneRigidbody.inertiaTensor = new Vector3(876f,947f,76f);
            PlaneRigidbody.inertiaTensorRotation = Quaternion.AngleAxis(-4.833f, Vector3.forward);
            // Specification At Cruise without Ground Effect
            Airspeed0 = 9.700f; // Magnitude of ground speed [m/s]
            alpha0 = 1.682f; // Angle of attack [deg]
            CDp0 = 0.018f; // Parasitic drag [-]
            Cmw0 = -0.164f; // Pitching momentum [-]
            CLMAX = 1.700f;
            // Wing
            Sw = 18.042f; // Wing area of wing [m^2]
            bw = 25.133f; // Wing span [m]
            cMAC = 0.757f; // Mean aerodynamic chord [m]
            aw = 0.108f; // Wing Lift Slope [1/deg]
            hw = (0.323f-0.250f); // Length between Wing a.c. and c.g.
            ew = 0.949f; // Wing efficiency
            AR = (bw*bw)/Sw; // Aspect Ratio
            // Tail
            Downwash = true; // Conventional Tail: True, T-Tail: False
            St = 1.375f; // Wing area of tail
            at = 0.076f; // Tail Lift Slope [1/deg]
            lt = 4.200f; // Length between Tail a.c. and c.g.
            deMAX = 10.000f; // Maximum elevator angle
            tau = 1.000f; // Control surface angle of attack effectiveness [-]
            VH = (St*lt)/(Sw*cMAC); // Tail Volume
            // Fin
            drMAX = 10.000f; // Maximum rudder angle            
            // Ground Effect
            CGEMIN = 0.215f; // Minimum Ground Effect Coefficient [-]
            // Stability derivatives
            Cyb = -0.002410f; // [1/deg]
            Cyp = -0.228437f; // [1/rad]
            Cyr = 0.090542f; // [1/rad]
            Cydr = 0.001908f; // [1/deg]
            Clb = -0.002002f; // [1/deg]
            Clp = -0.877559f; // [1/rad]
            Clr = 0.237651f; // [1/rad]
            Cldr = 0.000052f; // [1/deg]
            Cnb = -0.000059f; // [1/deg]
            Cnp = -0.142441f; // [1/rad]
            Cnr = -0.000491f; // [1/rad]
            Cndr = -0.000262f; // [1/deg]
        }else if(MyGameManeger.instance.PlaneName == "QX-19"){
            // Plane
            PlaneRigidbody.mass = 96.631f;
            PlaneRigidbody.centerOfMass = new Vector3(0f,0.294f,0f);
            PlaneRigidbody.inertiaTensor = new Vector3(991f,1032f,60f);
            PlaneRigidbody.inertiaTensorRotation = Quaternion.AngleAxis(-9.134f, Vector3.forward);
            // Specification At Cruise without Ground Effect
            Airspeed0 = 8.800f; // Magnitude of ground speed [m/s]
            alpha0 = 1.554f; // Angle of attack [deg]
            CDp0 = 0.019f; // Parasitic drag [-]
            Cmw0 = -0.170f; // Pitching momentum [-]
            CLMAX = 1.700f;
            // Wing
            Sw = 18.275f; // Wing area of wing [m^2]
            bw = 26.418f; // Wing span [m]
            cMAC = 0.736f; // Mean aerodynamic chord [m]
            aw = 0.105f; // Wing Lift Slope [1/deg]
            hw = (0.323f-0.250f); // Length between Wing a.c. and c.g.
            ew = 1.010f; // Wing efficiency
            AR = (bw*bw)/Sw; // Aspect Ratio
            // Tail
            Downwash = true; // Conventional Tail: True, T-Tail: False
            St = 1.548f; // Wing area of tail
            at = 0.082f; // Tail Lift Slope [1/deg]
            lt = 3.200f; // Length between Tail a.c. and c.g.
            deMAX = 10.000f; // Maximum elevator angle
            tau = 1.000f; // Control surface angle of attack effectiveness [-]
            VH = (St*lt)/(Sw*cMAC); // Tail Volume
            // Fin
            drMAX = 10.000f; // Maximum rudder angle            
            // Ground Effect
            CGEMIN = 0.361f; // Minimum Ground Effect Coefficient [-]
            // Stability derivatives
            Cyb = -0.005300f; // [1/deg]
            Cyp = -0.567798f; // [1/rad]
            Cyr = 0.225280f; // [1/rad]
            Cydr = 0.001721f; // [1/deg]
            Clb = -0.005118f; // [1/deg]
            Clp = -0.827488f; // [1/rad]
            Clr = 0.296796f; // [1/rad]
            Cldr = 0.000050f; // [1/deg]
            Cnb = -0.000808f; // [1/deg]
            Cnp = -0.165533f; // [1/rad]
            Cnr = 0.001675f; // [1/rad]
            Cndr = -0.000208f; // [1/deg]
        }else if(MyGameManeger.instance.PlaneName == "QX-20"){
            // Plane
            PlaneRigidbody.mass = 98.797f;
            PlaneRigidbody.centerOfMass = new Vector3(0f,0.29f,0f);
            PlaneRigidbody.inertiaTensor = new Vector3(1003f,1045f,58f);
            PlaneRigidbody.inertiaTensorRotation = Quaternion.AngleAxis(-9.112f, Vector3.forward);
            // Specification At Cruise without Ground Effect
            Airspeed0 = 9.600f; // Magnitude of ground speed [m/s]
            alpha0 = 1.459f; // Angle of attack [deg]
            CDp0 = 0.016f; // Parasitic drag [-]
            Cmw0 =-0.114f; // Pitching momentum [-]
            CLMAX = 1.700f;
            // Wing
            Sw = 18.816f; // Wing area of wing [m^2]
            bw = 26.679f; // Wing span [m]
            cMAC = 0.755f; // Mean aerodynamic chord [m]
            aw = 0.108f; // Wing Lift Slope [1/deg]
            hw = (0.323f-0.250f); // Length between Wing a.c. and c.g.
            ew = 0.986f; // Wing efficiency
            AR = (bw*bw)/Sw; // Aspect Ratio
            // Tail
            Downwash = false; // Conventional Tail: True, T-Tail: False
            St = 1.526f; // Wing area of tail
            at = 0.088f; // Tail Lift Slope [1/deg]
            lt = 3.200f; // Length between Tail a.c. and c.g.
            deMAX = 10.000f; // Maximum elevator angle
            tau = 1.000f; // Control surface angle of attack effectiveness [-]
            VH = (St*lt)/(Sw*cMAC); // Tail Volume
            // Fin
            drMAX = 10.000f; // Maximum rudder angle            
            // Ground Effect
            CGEMIN = 0.293f; // Minimum Ground Effect Coefficient [-]
            // Stability derivatives
            Cyb = -0.003555f; // [1/deg]
            Cyp = -0.455493f; // [1/rad]
            Cyr = 0.143466f; // [1/rad]
            Cydr = 0.000888f; // [1/deg]
            Clb = -0.004049f; // [1/deg]
            Clp = -0.829690f; // [1/rad]
            Clr = 0.227736f; // [1/rad]
            Cldr = 0.000016f; // [1/deg]
            Cnb = -0.000500f; // [1/deg]
            Cnp = -0.132307f; // [1/rad]
            Cnr = 0.000942f; // [1/rad]
            Cndr = -0.000106f; // [1/deg]
        }
    }

}
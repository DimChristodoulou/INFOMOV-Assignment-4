// RNG - Marsaglia's xor32
uint RandomUInt(uint seed)
{
	seed ^= seed << 13;
	seed ^= seed >> 17;
	seed ^= seed << 5;
	return seed;
}
float RandomFloat(uint seed) { return seed * 2.3283064365387e-10f; }

__kernel void dust( __global float4* posDirBuffer, 
                    __global float* peaksBuffer, 
                    __global uint* frameBuffer,
                    __global uint* frameChangeBuffer,
                    uint mapWidth,
                    uint mapHeight,
                    int sizeofPeaksBuffer,
                    uint frame){

    int id1d = get_global_id(0);
    uint seed = 0x12345678 + id1d + frame;
    float4 posDir = posDirBuffer[id1d];
   

    posDir.x += posDir.z;
    posDir.y += posDir.w;
    posDir.w *= 0.95f;
   
   

    if(posDir.x < 0){
        posDir.x = (float)(mapWidth - 1);
        seed = RandomUInt(seed);
        posDir.y = (float)(seed %mapHeight);
        seed = RandomUInt(seed);
        posDir.z = -1 - RandomFloat(seed) * 2;
        posDir.w = 0;
    }

   
    

    for (uint i = 0; i < sizeofPeaksBuffer; i++) {
   
        float topeakx = peaksBuffer[3 * i] - posDir.x;
        float topeaky = peaksBuffer[3 * i+1] - posDir.y;
        float2 topeak = (float2)(topeakx, topeaky);
        float g = peaksBuffer[3  * i + 2] * 0.02f / sqrt( dot( topeak, topeak ) );
        topeak = normalize( topeak );
        posDir.w -= topeak.y * g;
    }

   /* __local float _sharedMemory[44];
    int localId = get_local_id(0);
    float topeakx = peaksBuffer[3 * localId] - posx;
    float topeaky = peaksBuffer[3 * localId + 1] - posy;
    float2 topeak = (float2)(topeakx, topeaky);
    float g = peaksBuffer[3 * localId + 2] * 0.02f / sqrt(dot(topeak, topeak));
    topeak = normalize(topeak);
    _sharedMemory[localId] -= topeak.y * g;

    barrier(CLK_LOCAL_MEM_FENCE);

    if (localId = 0)
        for (int i = 0; i < 44; i++)
            diry += _sharedMemory[i];*/

   

    uint tempFrame = frameBuffer[id1d];
    uint tempFrameChange = frameChangeBuffer[id1d];
   
    seed = RandomUInt(seed);
    posDir.w += RandomFloat(seed) * 0.05f - 0.025f;
    posDirBuffer[id1d] = posDir;
    frameBuffer[id1d] = (tempFrame + tempFrameChange + 256) & 255;
   

}


